#!/usr/bin/env python3
"""Core end-to-end checks using the small dataset in tests/data."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import h5py

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BINARY = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))


def run(*args, expect_success=True):
    result = subprocess.run([str(arg) for arg in args], capture_output=True, text=True, timeout=90)
    if (result.returncode == 0) != expect_success:
        raise AssertionError(f"command: {args}\nstdout: {result.stdout}\nstderr: {result.stderr}")
    return result


def pair_rows(path):
    return sorted(line for line in path.read_text().splitlines() if not line.startswith('#'))


def bwa_mem2():
    configured = os.environ.get('RUST_HIC_BWA_MEM2')
    candidate = Path(configured) if configured else Path(shutil.which('bwa-mem2') or ROOT / 'envs/runhic/bin/bwa-mem2')
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise AssertionError('Install bwa-mem2 or set RUST_HIC_BWA_MEM2 to its executable')
    return candidate.resolve()


class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not BINARY.is_file():
            raise RuntimeError(f'Build rust-hic first: {BINARY} not found')
        cls.temporary = tempfile.TemporaryDirectory(prefix='rust-hic-tests-')
        cls.work = Path(cls.temporary.name)
        cls.index = cls.work / 'fastq-index'
        run(BINARY, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz', '--chunks', '3', '-o', cls.index)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_01_fastq_index(self):
        """Chunk boundaries cover each original mate exactly once, with matching IDs."""
        manifest = json.loads((self.index / 'manifest.json').read_text())
        self.assertEqual(manifest['chunks'], 3)
        self.assertTrue((self.index / 'r1.gzidx').is_file())
        self.assertTrue((self.index / 'r2.gzidx').is_file())
        mates = [gzip.decompress((DATA / f'R{i}.fastq.gz').read_bytes()) for i in (1, 2)]
        chunks = [[], []]
        for mate in (0, 1):
            key = f'r{mate + 1}_uncompressed_offset'
            offsets = [boundary[key] for boundary in manifest['boundaries']]
            self.assertEqual((offsets[0], offsets[-1]), (0, len(mates[mate])))
            for start, end in zip(offsets, offsets[1:]):
                lines = mates[mate][start:end].splitlines()
                self.assertGreater(len(lines), 0)
                self.assertEqual(len(lines) % 4, 0)
                chunks[mate].append([lines[i].split()[0].removesuffix(f'/{mate + 1}'.encode()) for i in range(0, len(lines), 4)])
            self.assertEqual(sum(map(len, chunks[mate])), 1024)
        self.assertEqual(chunks[0], chunks[1])
        default_index = self.work / 'default-index'
        run(BINARY, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz', '-o', default_index)
        self.assertEqual(json.loads((default_index / 'manifest.json').read_text())['chunks'], 1)

    def test_01_concatenated_gzip_members(self):
        """Index and read FASTQ inputs composed of multiple gzip members."""
        inputs = []
        expected = []
        for mate in (1, 2):
            data = gzip.decompress((DATA / f'R{mate}.fastq.gz').read_bytes())
            lines = data.splitlines(keepends=True)
            split = (len(lines) // 8) * 4
            path = self.work / f'concatenated-r{mate}.fastq.gz'
            path.write_bytes(gzip.compress(b''.join(lines[:split])) +
                             gzip.compress(b''.join(lines[split:])))
            inputs.append(path)
            expected.append(data.decode())
        index = self.work / 'concatenated-index'
        run(BINARY, 'index', *inputs, '--chunks', '1', '-o', index)
        for mate in (1, 2):
            result = run(BINARY, 'index', '--internal', 'stream', index, '0', str(mate))
            self.assertEqual(result.stdout, expected[mate - 1])

    def test_02_mapping_matches_bwa_mem2(self):
        """Built-in bwa-mem3 produces the same pairs as BWA-MEM2 on this dataset."""
        bwa = bwa_mem2()
        reference = self.work / 'genome.fa'
        shutil.copyfile(DATA / 'genome.fa', reference)
        run(bwa, 'index', reference)
        common = [BINARY, 'map', self.index, '-r', reference, '-x', reference, '--jobs', '2', '--threads-per-job', '2']
        run(*common, '--aligner', 'bwa-mem2', '--aligner-executable', bwa, '-o', self.work / 'bwa.pairs')
        run(*common, '-o', self.work / 'builtin.pairs')
        expected = pair_rows(self.work / 'bwa.pairs')
        actual = pair_rows(self.work / 'builtin.pairs')
        self.assertEqual(actual, expected, 'BWA-MEM2 differs')
        self.assertEqual(len(actual), 1020)
        self.assertEqual(sum(line.startswith('#columns:') for line in
                             (self.work / 'builtin.pairs').read_text().splitlines()), 1)
        run(*common, '--chunk', '1', '-o', self.work / 'single.pairs')
        single = pair_rows(self.work / 'single.pairs')
        self.assertTrue(single)
        self.assertTrue(set(single).issubset(actual))
        self.assertFalse((self.work / 'single.chunk-0001.pairs').exists())

    def test_03_pairs_processing(self):
        """Sort, filter, and deduplicate across chunks, then construct mcool."""
        clean = self.work / 'clean.pairs'
        run(BINARY, 'clean', DATA / 'chunk1.pairs', DATA / 'chunk2.pairs',
            '--chromsizes', DATA / 'chrom.sizes', '--fragments', DATA / 'fragments.bed',
            '--assembly', 'test', '--threads', '2', '-o', clean)
        lines = clean.read_text().splitlines()
        self.assertIn('#sorted: chr1-chr2-pos1-pos2', lines)
        rows = [line.split('\t') for line in lines if not line.startswith('#')]
        self.assertEqual(len(rows), 2)
        self.assertEqual({(r[1], int(r[2]), r[3], int(r[4])) for r in rows},
                         {('chr1', 1, 'chr1', 11), ('chr1', 20, 'chr2', 30)})
        original = clean.read_bytes()
        run(BINARY, 'clean', DATA / 'chunk1.pairs', '--chromsizes', DATA / 'chrom.sizes',
            '-o', clean, expect_success=False)
        self.assertEqual(clean.read_bytes(), original)
        # Sorting annotated pairs must retain all fields and exchange the two
        # fragment triplets whenever the two ends are flipped.
        annotated = self.work / 'annotated.pairs'
        annotated.write_text(
            '## pairs format v1.0.0\n'
            '#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type '
            'rfrag1 rfrag_start1 rfrag_end1 rfrag2 rfrag_start2 rfrag_end2\n'
            'flipped\tchr2\t30\tchr1\t20\t+\t-\tUR\t2\t21\t40\t1\t11\t20\n'
        )
        preserved = self.work / 'preserved.pairs'
        run(BINARY, 'clean', annotated, clean, '--chromsizes', DATA / 'chrom.sizes',
            '--assembly', 'test', '--threads', '2', '-o', preserved)
        rows14 = [line.split('\t') for line in pair_rows(preserved)]
        self.assertTrue(all(len(row) == 14 for row in rows14))
        flipped = next(row for row in rows14 if row[0] == 'flipped')
        self.assertEqual(flipped, ['flipped', 'chr1', '20', 'chr2', '30', '-', '+', 'RU',
                                   '1', '11', '20', '2', '21', '40'])
        # Reject unknown columns, rows that disagree with the header, and mixed
        # layouts across inputs; never silently drop an unrecognized field.
        invalid = self.work / 'invalid.pairs'
        for header, row in [
            ('readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type extra',
             'x\tchr1\t1\tchr1\t11\t+\t-\tUU\tvalue'),
            ('readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type',
             'x\tchr1\t1\tchr1\t11\t+\t-\tUU\tvalue'),
        ]:
            invalid.write_text(f'#columns: {header}\n{row}\n')
            run(BINARY, 'clean', invalid, '--chromsizes', DATA / 'chrom.sizes',
                '-o', self.work / 'invalid-clean.pairs', expect_success=False)
            run(BINARY, 'matrix', invalid, '--chromsizes', DATA / 'chrom.sizes',
                '-o', self.work / 'invalid.mcool', expect_success=False)
        run(BINARY, 'clean', DATA / 'chunk1.pairs', annotated,
            '--chromsizes', DATA / 'chrom.sizes', '-o', self.work / 'mixed.pairs',
            expect_success=False)
        run(BINARY, 'clean', annotated, '--chromsizes', DATA / 'chrom.sizes',
            '--fragments', DATA / 'fragments.bed', '-o', self.work / 'refilter.pairs',
            expect_success=False)
        # Exercise the separate merge path used by input_sorted: true.
        config = self.work / 'sorted-config.json'
        result_dir = self.work / 'sorted-run'
        split_lines = preserved.read_text().splitlines(keepends=True)
        header = ''.join(line for line in split_lines if line.startswith('#'))
        bodies = [line for line in split_lines if not line.startswith('#')]
        sorted_inputs = [self.work / f'sorted-{i}.pairs' for i in range(2)]
        for i, path in enumerate(sorted_inputs):
            path.write_text(header + ''.join(bodies[i::2]))
        config.write_text(json.dumps({
            'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
            'input_pairs': [str(path) for path in sorted_inputs], 'input_sorted': True,
            'output_dir': str(result_dir), 'scratch_dir': str(self.work),
            'resolutions': [5000], 'balance_max_iterations': 10,
        }))
        run(BINARY, 'run', config)
        self.assertEqual(pair_rows(result_dir / 'filtered-nodups.pairs'), pair_rows(preserved))
        mcool = self.work / 'contacts.mcool'
        run(BINARY, 'matrix', clean, '--chromsizes', DATA / 'chrom.sizes',
            '--assembly', 'test', '--threads', '2', '-o', mcool)
        with h5py.File(mcool, 'r') as f:
            self.assertEqual(len(f['resolutions']), 10)
            for resolution in ('5000', '10000', '5000000'):
                self.assertEqual(int(f[f'resolutions/{resolution}/pixels/count'][:].sum()), 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
