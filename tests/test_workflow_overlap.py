#!/usr/bin/env python3
"""Check overlapped mapping, filtering, global deduplication and failure handling."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

import h5py
import numpy as np

import layout_helpers as layout

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))


def run(*args):
    return subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=120)


def fake_chromap(path, bad=False):
    path.write_text('''#!/usr/bin/env python3
import concurrent.futures, pathlib, sys, time
args = sys.argv
r1 = pathlib.Path(args[args.index('-1') + 1])
r2 = pathlib.Path(args[args.index('-2') + 1])
chunk = int(r1.parent.name.split('-')[-1])
if BAD and chunk == 0:
    time.sleep(60)  # Never open the FIFOs; cancellation must release the producer.
def read(path):
    with open(path, 'rb') as f:
        return f.read()
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    a, b = list(pool.map(read, [r1, r2]))
assert a.count(b'\\n') == b.count(b'\\n')
if chunk == 0:
    time.sleep(3)
chrom = 'bad-chrom' if BAD and chunk == 1 else 'chr1'
output = pathlib.Path(args[args.index('-o') + 1])
header = ('## pairs format v1.0.0\\n'
          '#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\\n')
same = f'same{chunk}\\tchr1\\t{100 + chunk}\\tchr1\\t{200 + chunk}\\t+\\t-\\tUU\\n'
body = (f'dup{chunk}\\t{chrom}\\t2\\tchr2\\t20\\t+\\t-\\tUU\\n'
        f'unique{chunk}\\tchr1\\t{50 + 10 * chunk}\\tchr2\\t{100 + 10 * chunk}\\t+\\t-\\tUU\\n'
        + same)
# Chunk 0 becomes empty after restriction filtering.
output.write_text(header + (same if chunk == 0 else body))
print(f"Number of reads: {a.count(bytes([10])) // 2}.", file=sys.stderr)
'''.replace('BAD', repr(bad)))
    path.chmod(0o755)


class WorkflowOverlap(unittest.TestCase):
    def test_builtin_bwa_mem3_matches_sequential_mapping_and_clean(self):
        bwa = Path(os.environ.get('RUST_HIC_BWA_MEM2',
                                  shutil.which('bwa-mem2') or ROOT / 'envs/runhic/bin/bwa-mem2'))
        if not bwa.is_file():
            self.skipTest('bwa-mem2 is required to build the reference index')
        with tempfile.TemporaryDirectory(prefix='workflow-overlap-bwa-') as d:
            tmp = Path(d)
            reference = tmp / 'genome.fa'
            shutil.copyfile(DATA / 'genome.fa', reference)
            result = run(bwa, 'index', reference)
            self.assertEqual(result.returncode, 0, result.stderr)
            index = tmp / 'fastq-index'
            result = run(BIN, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz',
                         '--chunks', '3', '-o', index)
            self.assertEqual(result.returncode, 0, result.stderr)
            pairs = []
            for i in range(3):
                path = tmp / f'mapped-{i:04}.pairs'
                result = run(BIN, 'map', index, '-r', reference, '-x', reference,
                             '--aligner', 'bwa-mem3', '--jobs', '2', '--threads-per-job', '2',
                             '--chunk', str(i), '-o', path)
                self.assertEqual(result.returncode, 0, result.stderr)
                pairs.append(str(path))
            common = {
                'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
                'keep_pairs': True,
                'scratch_dir': str(tmp), 'resolutions': [1000], 'balance_max_iterations': 1,
                'resources': {'sort_jobs': 2, 'sort_threads': 1, 'sort_run_mib': 1,
                              'merge_workers': 2, 'dedup_workers': 2, 'matrix_threads': 2},
            }
            fastq = {
                'read1': str(DATA / 'R1.fastq.gz'), 'read2': str(DATA / 'R2.fastq.gz'),
                'reference': str(reference), 'aligner_index': str(reference),
                'index_dir': str(index), 'chunks': 3,
                'jobs': 2, 'threads_per_job': 2,
            }
            for name, inputs in [('sequential', {'input_pairs': pairs}), ('overlap', {'fastq': fastq})]:
                config = tmp / f'{name}.json'
                config.write_text(json.dumps({**common, **inputs, 'output_dir': str(tmp / name)}))
                result = run(BIN, 'run', config)
                self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('aligner: bwa-mem3',
                          (tmp / 'overlap' / 'config.resolved.yaml').read_text())
            self.assertEqual(layout.pairs(tmp / 'sequential').read_bytes(),
                             layout.pairs(tmp / 'overlap').read_bytes())
            with h5py.File(layout.mcool(tmp / 'overlap'), 'r') as actual, \
                 h5py.File(layout.mcool(tmp / 'sequential'), 'r') as expected:
                for column in ('bin1_id', 'bin2_id', 'count'):
                    np.testing.assert_array_equal(actual[f'resolutions/1000/pixels/{column}'][:],
                                                  expected[f'resolutions/1000/pixels/{column}'][:])

    def test_overlapped_chunks_match_sequential_clean_and_mcool(self):
        with tempfile.TemporaryDirectory(prefix='workflow-overlap-') as d:
            tmp = Path(d)
            aligner = tmp / 'fake-chromap'
            fake_chromap(aligner)
            reference_index = tmp / 'chromap.index'
            reference_index.touch()
            index = tmp / 'fastq-index'
            result = run(BIN, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz',
                         '--chunks', '3', '-o', index)
            self.assertEqual(result.returncode, 0, result.stderr)
            map_args = [BIN, 'map', index, '-r', DATA / 'genome.fa', '-x', reference_index,
                        '--aligner', 'chromap', '--aligner-executable', aligner,
                        '--jobs', '2', '--threads-per-job', '1']
            input_pairs = []
            for i in range(3):
                path = tmp / f'reference-mapped-{i:04}.pairs'
                result = run(*map_args, '--chunk', str(i), '-o', path)
                self.assertEqual(result.returncode, 0, result.stderr)
                input_pairs.append(str(path))
            combined = tmp / 'reference-mapped.pairs'
            result = run(*map_args, '-o', combined)
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = combined.read_text().splitlines()
            self.assertEqual(sum(line.startswith('#columns:') for line in lines), 1)
            rows = [line for line in lines if not line.startswith('#')]
            expected_rows = [line for path in input_pairs for line in Path(path).read_text().splitlines()
                             if not line.startswith('#')]
            self.assertEqual(rows, expected_rows)
            previous = combined.read_bytes()
            result = run(*map_args, '-o', combined)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(combined.read_bytes(), previous)
            common = {
                'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
                'keep_pairs': True,
                'scratch_dir': str(tmp), 'resolutions': [1000], 'balance_max_iterations': 1,
                'resources': {'sort_jobs': 2, 'sort_threads': 1, 'sort_run_mib': 1,
                              'merge_workers': 2, 'dedup_workers': 2, 'matrix_threads': 2},
            }
            fastq = {
                'read1': str(DATA / 'R1.fastq.gz'), 'read2': str(DATA / 'R2.fastq.gz'),
                'reference': str(DATA / 'genome.fa'), 'aligner_index': str(reference_index),
                'aligner': 'chromap', 'aligner_executable': str(aligner), 'index_dir': str(index),
                'chunks': 3, 'jobs': 2, 'threads_per_job': 1,
            }
            for use_fragments in (False, True):
                with self.subTest(fragments=use_fragments):
                    extra = {'fragments': str(DATA / 'fragments.bed')} if use_fragments else {}
                    sequential_dir = tmp / f'sequential-{use_fragments}'
                    overlap_dir = tmp / f'overlap-{use_fragments}'
                    sequential = {**common, **extra, 'input_pairs': input_pairs,
                                  'output_dir': str(sequential_dir)}
                    overlap = {**common, **extra, 'fastq': fastq, 'output_dir': str(overlap_dir)}
                    for label, config in [('sequential', sequential), ('overlap', overlap)]:
                        path = tmp / f'{label}-{use_fragments}.json'
                        path.write_text(json.dumps(config))
                        result = run(BIN, 'run', path)
                        self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(layout.pairs(overlap_dir).read_bytes(),
                                     layout.pairs(sequential_dir).read_bytes())
                    report = json.loads(layout.qc(overlap_dir, 'report.json').read_text())
                    self.assertEqual(report['input']['aligner'], 'chromap')
                    self.assertEqual(report['counts']['input_read_pairs'], 1024)
                    self.assertEqual(report['counts']['pairs_before_restriction'], 7)
                    self.assertEqual(report['restriction']['mode'], 'bed' if use_fragments else 'disabled')
                    self.assertTrue(all(chunk['selected_pairs'] is None for chunk in report['chunks']))
                    if use_fragments:
                        self.assertEqual(json.loads(layout.qc(overlap_dir, 'filter-stats.json').read_text()),
                                         json.loads(layout.qc(sequential_dir, 'filter-stats.json').read_text()))
                        stats = json.loads(layout.qc(overlap_dir, 'filter-stats.json').read_text())
                        self.assertEqual((stats['input_records'], stats['same_fragment']), (7, 3))
                        rows = [line for line in layout.pairs(overlap_dir).read_text().splitlines()
                                if not line.startswith('#')]
                        self.assertEqual(len(rows), 3)  # The duplicate spans chunks 1 and 2.
                    with h5py.File(layout.mcool(overlap_dir), 'r') as actual, \
                         h5py.File(layout.mcool(sequential_dir), 'r') as expected:
                        for column in ('bin1_id', 'bin2_id', 'count'):
                            self.assertEqual(actual[f'resolutions/1000/pixels/{column}'][:].tolist(),
                                             expected[f'resolutions/1000/pixels/{column}'][:].tolist())
                        for column in ('chrom', 'start', 'end'):
                            np.testing.assert_array_equal(actual[f'resolutions/1000/bins/{column}'][:],
                                                          expected[f'resolutions/1000/bins/{column}'][:])
                        np.testing.assert_allclose(actual['resolutions/1000/bins/weight'][:],
                                                   expected['resolutions/1000/bins/weight'][:],
                                                   equal_nan=True)
                    events = [json.loads(line) for line in layout.event_log(overlap_dir).read_text().splitlines()]
                    alignment_end = next(e['elapsed_seconds'] for e in events if e['event'] == 'alignment_done')
                    early = [e for e in events if e['event'] == 'sort_chunk_start'
                             and e['elapsed_seconds'] < alignment_end]
                    self.assertTrue(early, events)
                    if use_fragments:
                        self.assertTrue(any(e['event'] == 'filter_chunk_done'
                                            and e['elapsed_seconds'] < alignment_end for e in events), events)

            # A malformed completed chunk must fail the whole run, even while
            # another chunk is still mapping, and must not publish an mcool.
            bad_aligner = tmp / 'bad-chromap'
            fake_chromap(bad_aligner, bad=True)
            failed = {**common, 'fastq': {**fastq, 'aligner_executable': str(bad_aligner)},
                      'output_dir': str(tmp / 'failed')}
            config = tmp / 'failed.json'
            config.write_text(json.dumps(failed))
            started = time.monotonic()
            result = run(BIN, 'run', config)
            self.assertLess(time.monotonic() - started, 15, result.stderr)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('sort chunk 1', result.stderr)
            self.assertTrue((tmp / 'failed' / 'failure.json').is_file())
            failed_work = Path((tmp / 'failed' / 'work-dir.txt').read_text().strip())
            self.assertTrue(failed_work.is_dir())
            self.assertFalse(layout.has_mcool(tmp / 'failed'))

            bad_bed = tmp / 'incomplete.bed'
            bad_bed.write_text('chr1\t0\t25000\n')
            failed_filter = {**common, 'fastq': fastq, 'fragments': str(bad_bed),
                             'output_dir': str(tmp / 'failed-filter')}
            config = tmp / 'failed-filter.json'
            config.write_text(json.dumps(failed_filter))
            result = run(BIN, 'run', config)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('filter chunk ', result.stderr)
            self.assertTrue((tmp / 'failed-filter' / 'failure.json').is_file())
            self.assertFalse(layout.has_mcool(tmp / 'failed-filter'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
