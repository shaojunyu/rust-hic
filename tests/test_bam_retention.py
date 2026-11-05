#!/usr/bin/env python3
"""End-to-end coverage of optional unsorted BAM retention."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))
BWA = Path(os.environ.get('RUST_HIC_BWA_MEM2', ROOT / 'envs/runhic/bin/bwa-mem2'))
CHROMAP = Path(os.environ.get('RUST_HIC_CHROMAP', ROOT / 'envs/runhic/bin/chromap'))
SAMTOOLS = Path(os.environ.get('RUST_HIC_SAMTOOLS', ROOT / 'envs/runhic/bin/samtools'))


def command(*args):
    result = subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=120)
    if result.returncode:
        raise AssertionError(f'{args}: {result.stderr}')
    return result.stdout


@unittest.skipUnless(BWA.is_file() and CHROMAP.is_file() and SAMTOOLS.is_file(), 'aligners and samtools required')
class BamRetention(unittest.TestCase):
    def test_all_aligners_keep_all_emitted_records_without_changing_pairs(self):
        with tempfile.TemporaryDirectory(prefix='rust-hic-bam-') as d:
            tmp = Path(d)
            reference = tmp / 'genome.fa'
            shutil.copyfile(DATA / 'genome.fa', reference)
            command(BWA, 'index', reference)
            chromap_index = tmp / 'chromap.index'
            command(CHROMAP, '-i', '-r', reference, '-o', chromap_index)
            index = tmp / 'fastq-index'
            command(BIN, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz', '--chunks', '3', '-o', index)
            wrapper = tmp / 'capture-chromap'
            wrapper.write_text(f'''#!{sys.executable}\nimport os, pathlib, subprocess, sys\nexe = {str(CHROMAP)!r}\nif '--SAM' not in sys.argv:\n    os.execv(exe, [exe] + sys.argv[1:])\nchunk = pathlib.Path(sys.argv[sys.argv.index('-1') + 1]).parent.name\ncapture = pathlib.Path({str(tmp)!r}) / (chunk + '.sam')\nwith capture.open('wb') as output:\n    result = subprocess.run([exe] + sys.argv[1:], stdout=output)\nif result.returncode:\n    sys.exit(result.returncode)\nwith capture.open('rb') as data:\n    while block := data.read(1048576):\n        sys.stdout.buffer.write(block)\n''')
            wrapper.chmod(0o755)
            common = {'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
                      'keep_pairs': True, 'scratch_dir': str(tmp), 'resolutions': [1000],
                      'balance_max_iterations': 1,
                      'resources': {'sort_jobs': 2, 'sort_threads': 1, 'sort_run_mib': 1,
                                    'merge_workers': 2, 'dedup_workers': 2, 'matrix_threads': 2}}
            for aligner in ('bwa-mem3', 'bwa-mem2', 'chromap'):
                with self.subTest(aligner=aligner):
                    fastq = {'read1': str(DATA / 'R1.fastq.gz'), 'read2': str(DATA / 'R2.fastq.gz'),
                             'reference': str(reference), 'index_dir': str(index), 'chunks': 3,
                             'jobs': 2, 'threads_per_job': 2, 'aligner': aligner,
                             'aligner_index': str(chromap_index if aligner == 'chromap' else reference)}
                    if aligner != 'bwa-mem3':
                        fastq['aligner_executable'] = str(wrapper if aligner == 'chromap' else BWA)
                    results = []
                    for keep in (False, True):
                        out = tmp / f'{aligner}-{keep}'
                        config = tmp / f'{aligner}-{keep}.json'
                        yaml_option = keep and aligner == 'bwa-mem2'
                        config.write_text(json.dumps({**common, 'fastq': fastq, 'output_dir': str(out),
                                                      'keep_bam': yaml_option}))
                        command(BIN, 'run', config, *(['--keep-bam'] if keep and not yaml_option else []))
                        self.assertEqual((out / 'alignments.bam').exists(), keep)
                        self.assertFalse((out / 'work-dir.txt').exists())
                        results.append(out)
                    self.assertEqual((results[0] / 'filtered-nodups.pairs').read_bytes(),
                                     (results[1] / 'filtered-nodups.pairs').read_bytes())
                    with h5py.File(results[0] / 'contacts.mcool') as a, h5py.File(results[1] / 'contacts.mcool') as b:
                        for col in ('bin1_id', 'bin2_id', 'count'):
                            np.testing.assert_array_equal(a[f'resolutions/1000/pixels/{col}'][:], b[f'resolutions/1000/pixels/{col}'][:])
                    bam = results[1] / 'alignments.bam'
                    command(SAMTOOLS, 'quickcheck', bam)
                    header = command(SAMTOOLS, 'view', '-H', bam)
                    self.assertIn('SO:unsorted', header)
                    records = command(SAMTOOLS, 'view', bam).splitlines()
                    if aligner != 'chromap':
                        self.assertGreaterEqual(len(records), 2048)
                    else:
                        expected = []
                        for sam in sorted(tmp.glob('align-*.sam')):
                            chunk = int(sam.stem.split('-')[-1])
                            rows = [line.split('\t') for line in sam.read_text().splitlines() if not line.startswith('@')]
                            # Concurrent Chromap jobs must not exchange their
                            # low-memory temporary files or SAM records.
                            bounds = [0, 342, 683, 1024]
                            self.assertTrue(all(bounds[chunk] <= int(row[0][4:8]) < bounds[chunk + 1]
                                                for row in rows))
                            expected.extend(rows)
                        self.assertGreater(len(expected), 1000)
                        def normalized(row):
                            return (tuple(row[:11]), tuple(sorted(row[11:])))
                        self.assertEqual(sorted(normalized(r) for r in expected),
                                         sorted(normalized(r.split('\t')) for r in records))
                    self.assertTrue(all(len(row.split('\t')) >= 11 for row in records))
                    self.assertEqual(json.loads((results[1] / 'report.json').read_text())['outputs']['bam'], 'alignments.bam')


if __name__ == '__main__':
    unittest.main(verbosity=2)
