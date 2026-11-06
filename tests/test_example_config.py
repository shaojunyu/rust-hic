#!/usr/bin/env python3
"""Exercise the embedded example configuration through the public CLI."""
import gzip
import json
from collections import Counter
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

    def test_minimal_example_is_embedded(self):
        result = command(BIN, 'example-config', '--minimal')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, '')
        self.assertEqual(result.stdout, (ROOT / 'docs/example.minimal.config.yaml').read_text())

    @unittest.skipUnless(BWA.is_file(), 'bwa-mem2 is needed to build a reference index')
    def test_minimal_fastq_example_runs_with_defaults(self):
        result = command(BWA, 'index', self.work / 'hg38.fa')
        self.assertEqual(result.returncode, 0, result.stderr)
        for mate in ('R1', 'R2'):
            shutil.copy(DATA / f'{mate}.fastq.gz', self.work / f'{mate}.fastq.gz')
        example = command(BIN, 'example-config', '--minimal')
        self.assertEqual(example.returncode, 0, example.stderr)
        config = self.work / 'minimal.yaml'
        config.write_text(example.stdout)
        check = command(BIN, 'run', config, '--check')
        self.assertEqual(check.returncode, 0, check.stderr)
        self.assertFalse((self.work / 'new-result').exists())
        result = command(BIN, 'run', config)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.work / 'new-result'
        self.assertTrue((out / 'contacts.mcool').is_file())
        self.assertTrue((out / 'report.html').is_file())
        self.assertFalse((out / 'filtered-nodups.pairs').exists())
        self.assertFalse((out / 'alignments.bam').exists())
        resolved = (out / 'config.resolved.yaml').read_text()
        self.assertIn('aligner: bwa-mem3', resolved)
        self.assertIn('index_threads_per_mate: 8', resolved)
        self.assertIn('balance_max_iterations: 1000', resolved)
        self.assertIn('fragments: null', resolved)

    def run_example(self, text, mismatch=3):
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
        text = text.replace('keep_pairs: false', 'keep_pairs: true')
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
                     'filtered-nodups.pairs.gz.px2', 'status.json', 'config.resolved.yaml',
                     'report.html', 'report.json'):
            self.assertTrue((output / name).is_file(), name)
        resolved = (output / 'config.resolved.yaml').read_text()
        self.assertIn('compressed_pairs: true', resolved)
        self.assertIn('balance_max_iterations: 2', resolved)
        self.assertIn(f'dedup_max_mismatch: {mismatch}', resolved)
        self.assertIn('sort_run_mib: 1', resolved)
        self.assertEqual(json.loads((output / 'status.json').read_text())['status'], 'complete')
        report = json.loads((output / 'report.json').read_text())
        counts = report['counts']
        self.assertEqual(counts['pairs_before_dedup'] - counts['duplicates'], counts['final_contacts'])
        self.assertEqual(counts['cis'] + counts['trans'], counts['final_contacts'])
        self.assertEqual([item['bp'] for item in report['resolutions']], [1000, 2000])
        for item in report['resolutions']:
            if item['usable_bins'] == 0:
                self.assertIsNone(item['converged'])
        self.assertIn('Resolution &amp; balance', (output / 'report.html').read_text())
        records = [line.split('\t') for line in (output / 'filtered-nodups.pairs').read_text().splitlines()
                   if not line.startswith('#')]
        self.assertEqual(counts['pair_types'], dict(Counter(row[7] for row in records)))
        self.assertEqual(sum(counts['strand_orientations'].values()), len(records))
        self.assertEqual(sum(item['count'] for item in counts['cis_distance_bins']), counts['cis'])
        return resolved

    def test_example_with_pairs_and_digest(self):
        shutil.copy(DATA / 'chunk1.pairs', self.work / 'sample1.pairs')
        shutil.copy(DATA / 'chunk2.pairs', self.work / 'sample2.pairs')
        text = self.example.replace('# fragments: hg38.Arima.bed', 'digest:\n  enzyme: Arima\n  reference: hg38.fa')
        start, end = text.index('fastq:\n'), text.index('# To start from uncompressed pairs')
        text = text[:start] + 'input_pairs: [sample1.pairs, sample2.pairs]\n\n' + text[end:]
        resolved = self.run_example(text)
        self.assertIn('enzyme: Arima', resolved)
        self.assertIn('sample1.pairs', resolved)
        report = json.loads((self.work / 'new-result/report.json').read_text())
        self.assertEqual(report['restriction']['mode'], 'arima_digest')
        self.assertEqual(report['restriction']['stats']['input_records'],
                         report['counts']['pairs_before_restriction'])
        self.assertIsNone(report['counts']['input_read_pairs'])

    def test_example_without_restriction_filtering(self):
        shutil.copy(DATA / 'chunk1.pairs', self.work / 'sample1.pairs')
        start = self.example.index('fastq:\n')
        end = self.example.index('# To start from uncompressed pairs')
        text = self.example[:start] + 'input_pairs: [sample1.pairs]\n\n' + self.example[end:]
        resolved = self.run_example(text)
        self.assertIn('fragments: null', resolved)
        self.assertIn('digest: null', resolved)
        self.assertFalse((self.work / 'new-result' / 'filter-stats.json').exists())
        pairs = (self.work / 'new-result' / 'filtered-nodups.pairs').read_text()
        self.assertIn('#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n', pairs)
        report = json.loads((self.work / 'new-result/report.json').read_text())
        self.assertEqual(report['restriction']['mode'], 'disabled')
        self.assertIsNone(report['restriction']['stats'])
        dest = self.work / 'regenerated'
        dest.mkdir()
        result = command(BIN, 'run', 'report', self.work / 'new-result', '--output', dest)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((dest / 'report.json').read_text())['counts'], report['counts'])
        result = command(BIN, 'run', 'report', self.work / 'new-result', '--output', dest)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('already exists', result.stderr)
        # Older run outputs lack the detailed distributions; never invent them.
        stats_path = self.work / 'new-result/contact-stats.json'
        stats = json.loads(stats_path.read_text())
        for name in ('pair_types', 'strand_orientations', 'cis_distance_bins'):
            del stats[name]
        stats_path.write_text(json.dumps(stats))
        old_dest = self.work / 'old-report'
        old_dest.mkdir()
        result = command(BIN, 'run', 'report', self.work / 'new-result', '--output', old_dest)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIsNone(json.loads((old_dest / 'report.json').read_text())['counts']['cis_distance_bins'])
        self.assertIn('were not recorded', (old_dest / 'report.html').read_text())
        recalc_dest = self.work / 'recalculated'
        recalc_dest.mkdir()
        result = command(BIN, 'run', 'report', self.work / 'new-result', '--output', recalc_dest,
                         '--recompute-contact-stats')
        self.assertEqual(result.returncode, 0, result.stderr)
        recalculated = json.loads((recalc_dest / 'report.json').read_text())
        self.assertEqual(recalculated['counts']['cis_distance_bins'], report['counts']['cis_distance_bins'])
        self.assertNotIn('cis_distance_bins', json.loads(stats_path.read_text()))

    def test_distance_boundaries_and_escaped_metadata(self):
        distances = [0, 999, 1000, 1999, 2000, 4999, 5000, 9999, 10000, 19999, 20000]
        header = ('## pairs format v1.0.0\n'
                  '#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n')
        rows = [f'read{i}\tchr1\t1000\tchr1\t{1000 + distance}\t+\t-\tUU\n'
                for i, distance in enumerate(distances)]
        (self.work / 'sample1.pairs').write_text(header + ''.join(rows))
        start = self.example.index('fastq:\n')
        end = self.example.index('# To start from uncompressed pairs')
        text = self.example[:start] + 'input_pairs: [sample1.pairs]\n\n' + self.example[end:]
        text = text.replace('assembly: hg38', "assembly: 'hg38<&'")
        text = text.replace('dedup_max_mismatch: 3', 'dedup_max_mismatch: 0')
        self.run_example(text, mismatch=0)
        report = json.loads((self.work / 'new-result/report.json').read_text())
        records = [line.split('\t') for line in (self.work / 'new-result/filtered-nodups.pairs').read_text().splitlines()
                   if not line.startswith('#')]
        actual_distances = [abs(int(row[2]) - int(row[4])) for row in records if row[1] == row[3]]
        self.assertCountEqual(actual_distances, distances)
        for bin in report['counts']['cis_distance_bins']:
            self.assertEqual(bin['count'], sum(bin['start_bp'] <= distance and
                                             (bin['end_bp'] is None or distance < bin['end_bp'])
                                             for distance in actual_distances))
        page = (self.work / 'new-result/report.html').read_text()
        self.assertIn('hg38&lt;&amp;', page)
        self.assertNotIn('hg38<&', page)

    def test_pairs_retention_defaults_off_and_cli_enables_it(self):
        base = {
            'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
            'input_pairs': [str(DATA / 'chunk1.pairs')], 'scratch_dir': str(self.work),
            'resolutions': [1000], 'balance_max_iterations': 1,
            'resources': {'sort_jobs': 1, 'sort_threads': 1, 'sort_run_mib': 1,
                          'merge_workers': 1, 'dedup_workers': 1, 'matrix_threads': 1},
        }
        for retain in (False, True):
            out = self.work / f'retention-{retain}'
            config = self.work / f'retention-{retain}.json'
            config.write_text(json.dumps({**base, 'output_dir': str(out)}))
            result = command(BIN, 'run', config, *(['--keep-pairs'] if retain else []))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((out / 'contacts.mcool').is_file())
            self.assertTrue((out / 'report.html').is_file())
            self.assertEqual((out / 'filtered-nodups.pairs').exists(), retain)
            self.assertFalse((out / 'work-dir.txt').exists())
            report = json.loads((out / 'report.json').read_text())
            self.assertEqual(report['outputs']['pairs'] is not None, retain)
            if not retain:
                target = self.work / 'recompute-without-pairs'
                target.mkdir()
                result = command(BIN, 'run', 'report', out, '--output', target, '--recompute-contact-stats')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('final pairs were not retained', result.stderr)
        invalid = self.work / 'pairs-bam.json'
        invalid.write_text(json.dumps({**base, 'output_dir': str(self.work / 'invalid-bam')}))
        result = command(BIN, 'run', invalid, '--check', '--keep-bam')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires FASTQ input', result.stderr)
        self.assertFalse(list(self.work.glob('rust-hic-*')))

    @unittest.skipUnless(BWA.is_file(), 'bwa-mem2 is needed to build a reference index')
    def test_example_with_builtin_bwa_mem3(self):
        result = command(BWA, 'index', self.work / 'hg38.fa')
        self.assertEqual(result.returncode, 0, result.stderr)
        for mate in ('R1', 'R2'):
            shutil.copy(DATA / f'{mate}.fastq.gz', self.work / f'{mate}.fastq.gz')
        with gzip.open(self.work / 'R1.fastq.gz', 'rt') as reads:
            count = sum(1 for _ in reads) // 4
        text = self.example.replace('  # chunks: 16', '  chunks: 3')
        text = text.replace('# fragments: hg38.Arima.bed', 'fragments: hg38.Arima.bed')
        text = text.replace('  jobs: 4', '  jobs: 2')
        text = text.replace('  threads_per_job: 8', '  threads_per_job: 1')
        text = text.replace('  index_threads_per_mate: 8', '  index_threads_per_mate: 1')
        text = text.replace('  # expected_pairs: 1000000', f'  expected_pairs: {count}')
        resolved = self.run_example(text)
        self.assertIn('aligner: bwa-mem3', resolved)
        self.assertIn(f'expected_pairs: {count}', resolved)
        self.assertIn('index_threads_per_mate: 1', resolved)
        self.assertIn('chunks: 3', resolved)
        report = json.loads((self.work / 'new-result/report.json').read_text())
        self.assertEqual(report['restriction']['mode'], 'bed')
        self.assertEqual(report['counts']['input_read_pairs'], count)
        self.assertEqual(len(report['chunks']), 3)
        self.assertEqual(sum(item['selected_pairs'] for item in report['chunks']),
                         report['counts']['pairs_before_restriction'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
