#!/usr/bin/env python3
"""Public preflight behavior and actionable errors."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='rust-hic-preflight-')
        self.work = Path(self.tmp.name)
        self.output = self.work / 'result'
        self.config = {
            'assembly': 'test',
            'chromsizes': str(DATA / 'chrom.sizes'),
            'input_pairs': [str(DATA / 'chunk1.pairs')],
            'output_dir': str(self.output),
            'scratch_dir': str(self.work),
        }

    def tearDown(self):
        self.tmp.cleanup()

    def command(self, *args):
        return subprocess.run([str(BIN), *map(str, args)], capture_output=True, text=True, timeout=20)

    def check(self, check_only=True):
        path = self.work / 'config.json'
        path.write_text(json.dumps(self.config))
        return self.command('run', path, *(['--check'] if check_only else []))

    def test_valid_check_does_not_start_pipeline(self):
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Preflight passed', result.stdout)
        self.assertFalse(self.output.exists())
        self.assertFalse(list(self.work.glob('rust-hic-*')))

    def test_invalid_sizes_fail_before_output_is_created(self):
        sizes = self.work / 'bad.sizes'
        sizes.write_text('chr1 100\nchr1 200\n')
        self.config['chromsizes'] = str(sizes)
        for check_only in (True, False):
            result = self.check(check_only)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('bad.sizes:2', result.stderr)
            self.assertFalse(self.output.exists())

    def test_bad_bed_and_missing_pairs_fail_early(self):
        bed = self.work / 'bad.bed'
        bed.write_text('chr1 0 10\nchr1 12 20\n')
        self.config['fragments'] = str(bed)
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('bad.bed', result.stderr)
        self.assertIn('not contiguous', result.stderr)
        del self.config['fragments']
        self.config['input_pairs'] = [str(self.work / 'missing.pairs')]
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('missing.pairs', result.stderr)
        self.assertFalse(self.output.exists())

    def test_invalid_resource_names_are_reported(self):
        self.config['resources'] = {'merge_workers': 65}
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('resources.merge_workers', result.stderr)
        self.assertFalse(self.output.exists())

    def test_compressed_pairs_input_has_actionable_error(self):
        path = self.work / 'input.pairs.gz'
        path.write_bytes(gzip.compress((DATA / 'chunk1.pairs').read_bytes()))
        self.config['input_pairs'] = [str(path)]
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('pairs input must be uncompressed', result.stderr)
        self.assertIn('input.pairs.gz', result.stderr)
        self.assertFalse(self.output.exists())

    def test_sorted_headers_must_match_and_output_must_be_new(self):
        second = self.work / 'second.pairs'
        second.write_text((DATA / 'chunk2.pairs').read_text().replace('## pairs format', '#extra\n## pairs format', 1))
        self.config['input_pairs'].append(str(second))
        self.config['input_sorted'] = True
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('headers do not match', result.stderr)
        self.config['input_pairs'].pop()
        self.output.mkdir()
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('output already exists', result.stderr)

    def test_non_gzip_fastq_does_not_create_index_directory(self):
        bad = self.work / 'not-gzip.fastq.gz'
        bad.write_bytes(b'not gzip')
        output = self.work / 'index'
        result = self.command('index', bad, DATA / 'R2.fastq.gz', '-o', output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not gzip-compressed', result.stderr)
        self.assertIn('not-gzip.fastq.gz', result.stderr)
        self.assertFalse(output.exists())

    def test_hardlinked_mates_are_rejected_before_indexing(self):
        r1 = self.work / 'r1.gz'
        r2 = self.work / 'r2.gz'
        shutil.copyfile(DATA / 'R1.fastq.gz', r1)
        os.link(r1, r2)
        output = self.work / 'index'
        result = self.command('index', r1, r2, '-o', output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('same FASTQ file', result.stderr)
        self.assertFalse(output.exists())

    def test_reused_index_fingerprint_is_checked_before_output(self):
        index = self.work / 'index'
        result = self.command('index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz',
                              '--chunks', '1', '--threads-per-mate', '1', '-o', index)
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest_path = index / 'manifest.json'
        manifest = json.loads(manifest_path.read_text())
        manifest['inputs'][0]['edge_sha256'] = '0' * 64
        manifest_path.write_text(json.dumps(manifest))
        aligner_index = self.work / 'chromap.index'
        aligner_index.touch()
        del self.config['input_pairs']
        self.config['fastq'] = {
            'read1': str(DATA / 'R1.fastq.gz'), 'read2': str(DATA / 'R2.fastq.gz'),
            'reference': str(DATA / 'genome.fa'), 'aligner_index': str(aligner_index),
            'aligner': 'chromap', 'aligner_executable': '/bin/true', 'index_dir': str(index),
        }
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('input changed', result.stderr)
        self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
