#!/usr/bin/env python3
"""End-to-end checks for the paired FASTQ stream sent to Chromap."""
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))
CHROMAP = os.environ.get('RUST_HIC_CHROMAP') or shutil.which('chromap')


def run(*args):
    return subprocess.run([str(x) for x in args], capture_output=True, text=True, timeout=90)


@unittest.skipUnless(CHROMAP, 'set RUST_HIC_CHROMAP or install chromap')
class ChromapPairing(unittest.TestCase):
    def test_multiple_chunks_are_combined_in_order(self):
        with tempfile.TemporaryDirectory(prefix='chromap-map-combined-') as d:
            tmp = Path(d)
            reference_index = tmp / 'genome.index'
            result = run(CHROMAP, '-i', '-r', DATA / 'genome.fa', '-o', reference_index)
            self.assertEqual(result.returncode, 0, result.stderr)
            index = tmp / 'fastq-index'
            result = run(BIN, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz',
                         '--chunks', '3', '-o', index)
            self.assertEqual(result.returncode, 0, result.stderr)
            args = [BIN, 'map', index, '-r', DATA / 'genome.fa', '-x', reference_index,
                    '--aligner', 'chromap', '--aligner-executable', CHROMAP,
                    '--jobs', '2', '--threads-per-job', '1']
            combined = tmp / 'combined.pairs'
            result = run(*args, '-o', combined)
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = combined.read_text().splitlines()
            self.assertEqual(sum(line.startswith('#columns:') for line in lines), 1)
            expected = []
            for i in range(3):
                part = tmp / f'part-{i}.pairs'
                result = run(*args, '--chunk', str(i), '-o', part)
                self.assertEqual(result.returncode, 0, result.stderr)
                expected.extend(line for line in part.read_text().splitlines() if not line.startswith('#'))
            self.assertEqual([line for line in lines if not line.startswith('#')], expected)

    def test_valid_and_mismatched_internal_ids(self):
        with tempfile.TemporaryDirectory(prefix='chromap-pairing-') as tmp:
            tmp = Path(tmp)
            chromap_index = tmp / 'genome.index'
            result = run(CHROMAP, '-i', '-r', DATA / 'genome.fa', '-o', chromap_index)
            self.assertEqual(result.returncode, 0, result.stderr)
            r1 = DATA / 'R1.fastq.gz'
            original = gzip.decompress((DATA / 'R2.fastq.gz').read_bytes()).splitlines(keepends=True)
            self.assertGreater(len(original), 2000)
            # Keep first and last pair intact: indexing checks boundaries only.
            mismatched = original.copy()
            mismatched[4 * 500] = mismatched[4 * 500].replace(b'@', b'@mismatch_', 1)
            bad_r2 = tmp / 'bad-r2.fastq.gz'
            with gzip.open(bad_r2, 'wb') as f:
                f.write(b''.join(mismatched))
            missing_r2 = tmp / 'missing-r2.fastq.gz'
            with gzip.open(missing_r2, 'wb') as f:
                f.write(b''.join(original[:4 * 500] + original[4 * 501:]))
            for label, r2 in [('valid', DATA / 'R2.fastq.gz'), ('mismatch', bad_r2),
                              ('missing', missing_r2)]:
                index = tmp / (label + '-index')
                result = run(BIN, 'index', r1, r2, '--chunks', '1', '-o', index)
                self.assertEqual(result.returncode, 0, result.stderr)
                result = run(BIN, 'map', index, '-r', DATA / 'genome.fa', '-x', chromap_index,
                             '--aligner', 'chromap', '--aligner-executable', CHROMAP,
                             '--jobs', '1', '--threads-per-job', '2', '-o', tmp / (label + '.pairs'))
                if label == 'valid':
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertTrue((tmp / 'valid.pairs').is_file())
                else:
                    self.assertNotEqual(result.returncode, 0, result.stderr)
                    self.assertIn('FASTQ mate IDs differ', result.stderr)
            failing_aligner = tmp / 'failing-chromap'
            failing_aligner.write_text('#!/bin/sh\nexit 1\n')
            failing_aligner.chmod(0o755)
            result = run(BIN, 'map', tmp / 'valid-index', '-r', DATA / 'genome.fa',
                         '-x', chromap_index, '--aligner', 'chromap',
                         '--aligner-executable', failing_aligner, '--jobs', '1',
                         '--threads-per-job', '2', '-o', tmp / 'failed.pairs')
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((tmp / 'failed.pairs').exists())

            # A consumer can read an entire mate before the other. Both inputs
            # exceed the FIFO and in-memory buffer, so this exercises spilling.
            expected = [hashlib.sha256(gzip.decompress((DATA / f'R{i}.fastq.gz').read_bytes())).hexdigest()
                        for i in (1, 2)]
            for first in (1, 2):
                consumer = tmp / f'sequential-chromap-{first}'
                consumer.write_text(
                    '#!/usr/bin/env python3\n'
                    'import hashlib, sys\n'
                    'a=sys.argv\n'
                    'paths=[a[a.index("-1")+1], a[a.index("-2")+1]]\n'
                    'readers=[open(p,"rb",buffering=0) for p in paths]\n'
                    f'expected={expected!r}\n'
                    f'order={[first - 1, 2 - first]!r}\n'
                    'for i in order:\n'
                    '    h=hashlib.sha256()\n'
                    '    while True:\n'
                    '        data=readers[i].read(65536)\n'
                    '        if not data: break\n'
                    '        h.update(data)\n'
                    '    assert h.hexdigest() == expected[i]\n'
                    'print("Number of reads: 2048.", file=sys.stderr)\n'
                    'open(a[a.index("-o")+1],"w").write('
                    '"## pairs format v1.0.0\\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\\n")\n'
                )
                consumer.chmod(0o755)
                result = run(BIN, 'map', tmp / 'valid-index', '-r', DATA / 'genome.fa',
                             '-x', chromap_index, '--aligner', 'chromap',
                             '--aligner-executable', consumer, '--jobs', '1',
                             '--threads-per-job', '1', '-o', tmp / f'sequential-{first}.pairs')
                self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
