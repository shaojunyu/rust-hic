#!/usr/bin/env python3
"""Exercise the embedded example configuration through the public CLI."""
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
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic')).resolve()
BWA = Path(os.environ.get('RUST_HIC_BWA_MEM2',
                          shutil.which('bwa-mem2') or ROOT / 'envs/runhic/bin/bwa-mem2'))


def command(*args):
    return subprocess.run([str(arg) for arg in args], capture_output=True, text=True, timeout=120)


class ExampleConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='rust-hic-example-')
        self.work = Path(self.tmp.name)
        result = command(BIN, 'example-config')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, '')
        self.assertEqual(result.stdout, (ROOT / 'docs/example.config.yaml').read_text())
        self.example = result.stdout
        shutil.copy(DATA / 'chrom.sizes', self.work / 'hg38.chrom.sizes')
        shutil.copy(DATA / 'genome.fa', self.work / 'hg38.fa')
        shutil.copy(DATA / 'fragments.bed', self.work / 'hg38.Arima.bed')

    def tearDown(self):
        self.tmp.cleanup()

    def run_example(self, text):
        # Keep the tiny fixture fast and exercise nondefault settings.
        text = text.replace('  merge_workers: 6', '  merge_workers: 2')
        text = text.replace('  dedup_workers: 6', '  dedup_workers: 2')
        text = text.replace('  matrix_threads: 6', '  matrix_threads: 2')
        text = text.replace('  sort_jobs: 4', '  sort_jobs: 2')
        text = text.replace('  sort_threads: 2', '  sort_threads: 1')
        text = text.replace('  sort_run_mib: 256', '  sort_run_mib: 1')
        text = text.replace('resolutions: [5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000]',
                            'resolutions: [1000, 2000]')
        text = text.replace('balance_max_iterations: 1000', 'balance_max_iterations: 2')
        text = text.replace('compressed_pairs: false', 'compressed_pairs: true')
        config = self.work / 'config.yaml'
        config.write_text(text)
        result = command(BIN, 'run', config, '--check')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Preflight passed', result.stdout)
        self.assertFalse((self.work / 'new-result').exists())
        result = command(BIN, 'run', config)
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.work / 'new-result'
        for name in ('contacts.mcool', 'filtered-nodups.pairs', 'filtered-nodups.pairs.gz',
                     'filtered-nodups.pairs.gz.px2', 'status.json', 'config.resolved.yaml'):
            self.assertTrue((output / name).is_file(), name)
        resolved = (output / 'config.resolved.yaml').read_text()
        self.assertIn('compressed_pairs: true', resolved)
        self.assertIn('balance_max_iterations: 2', resolved)
        self.assertIn('dedup_max_mismatch: 3', resolved)
        self.assertIn('sort_run_mib: 1', resolved)
        self.assertEqual(json.loads((output / 'status.json').read_text())['status'], 'complete')
        return resolved

    def test_example_with_pairs_and_digest(self):
        shutil.copy(DATA / 'chunk1.pairs', self.work / 'sample1.pairs')
        shutil.copy(DATA / 'chunk2.pairs', self.work / 'sample2.pairs')
        text = self.example.replace('fragments: hg38.Arima.bed', 'digest:\n  enzyme: Arima\n  reference: hg38.fa')
        start, end = text.index('fastq:\n'), text.index('# To start from uncompressed pairs')
        text = text[:start] + 'input_pairs: [sample1.pairs, sample2.pairs]\n\n' + text[end:]
        resolved = self.run_example(text)
        self.assertIn('enzyme: Arima', resolved)
        self.assertIn('sample1.pairs', resolved)

    @unittest.skipUnless(BWA.is_file(), 'bwa-mem2 is needed to build a reference index')
    def test_example_with_builtin_bwa_mem3(self):
        result = command(BWA, 'index', self.work / 'hg38.fa')
        self.assertEqual(result.returncode, 0, result.stderr)
        for mate in ('R1', 'R2'):
            shutil.copy(DATA / f'{mate}.fastq.gz', self.work / f'{mate}.fastq.gz')
        with gzip.open(self.work / 'R1.fastq.gz', 'rt') as reads:
            count = sum(1 for _ in reads) // 4
        text = self.example.replace('  # chunks: 16', '  chunks: 3')
        text = text.replace('  jobs: 4', '  jobs: 2')
        text = text.replace('  threads_per_job: 8', '  threads_per_job: 1')
        text = text.replace('  index_threads_per_mate: 8', '  index_threads_per_mate: 1')
        text = text.replace('  # expected_pairs: 1000000', f'  expected_pairs: {count}')
        resolved = self.run_example(text)
        self.assertIn('aligner: bwa-mem3', resolved)
        self.assertIn(f'expected_pairs: {count}', resolved)
        self.assertIn('index_threads_per_mate: 1', resolved)
        self.assertIn('chunks: 3', resolved)


if __name__ == '__main__':
    unittest.main(verbosity=2)
