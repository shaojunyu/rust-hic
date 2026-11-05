#!/usr/bin/env python3
"""Exercise Slurm submission and workers with a local scheduler stub."""
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'tests' / 'data'
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic')).resolve()
BWA = Path(os.environ.get('RUST_HIC_BWA_MEM2', ROOT / 'envs/runhic/bin/bwa-mem2')).resolve()
SAMTOOLS = Path(os.environ.get('RUST_HIC_SAMTOOLS', ROOT / 'envs/runhic/bin/samtools')).resolve()
CHROMAP = Path(os.environ.get('RUST_HIC_CHROMAP', ROOT / 'envs/runhic/bin/chromap')).resolve()


def call(*args, env=None):
    return subprocess.run([str(a) for a in args], env=env, text=True, capture_output=True, timeout=180)


def write_scheduler(directory):
    sbatch = directory / 'sbatch'
    sbatch.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
args = sys.argv[1:]
if args == ['--version']:
    print('slurm 25.11.8')
    raise SystemExit(0)
script = pathlib.Path(args[-1]).name
if os.environ.get('FAIL_SBATCH') == script:
    print('deliberate scheduler failure', file=sys.stderr)
    raise SystemExit(1)
root = pathlib.Path(os.environ['FAKE_SLURM_DIR'])
counter = root / 'counter'
n = int(counter.read_text()) + 1 if counter.exists() else 1001
counter.write_text(str(n))
with (root / 'submissions.jsonl').open('a') as f:
    f.write(json.dumps({'id': str(n), 'args': args}) + '\\n')
print(n)
''')
    scontrol = directory / 'scontrol'
    scontrol.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
args = sys.argv[1:]
if args == ['--version']:
    if os.environ.get('FAIL_SCONTROL') == 'version':
        print('deliberate scontrol version failure', file=sys.stderr)
        raise SystemExit(1)
    print('slurm 25.11.8')
    raise SystemExit(0)
root = pathlib.Path(os.environ['FAKE_SLURM_DIR'])
with (root / 'control.jsonl').open('a') as f:
    f.write(json.dumps(args) + '\\n')
if os.environ.get('FAIL_SCONTROL') == args[0]:
    print('deliberate scontrol failure', file=sys.stderr)
    raise SystemExit(1)
if args[0] == 'update':
    job = next(x.split('=', 1)[1] for x in args if x.startswith('JobId='))
    name = next(x.split('=', 1)[1] for x in args if x.startswith('JobName='))
    path = root / 'names.json'
    names = json.loads(path.read_text()) if path.exists() else {}
    names[job] = name
    path.write_text(json.dumps(names))
elif args[0] == 'release':
    submissions = (root / 'submissions.jsonl').read_text().splitlines()
    if len(submissions) < 3:
        print('array released before final job submission', file=sys.stderr)
        raise SystemExit(1)
else:
    raise SystemExit(2)
''')
    squeue = directory / 'squeue'
    squeue.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
root = pathlib.Path(os.environ['FAKE_SLURM_DIR'])
p = root / 'active'
ids = p.read_text().splitlines() if p.exists() else []
path = root / 'names.json'
names = json.loads(path.read_text()) if path.exists() else {}
submissions = root / 'submissions.jsonl'
if submissions.exists():
    for line in submissions.read_text().splitlines():
        item = json.loads(line)
        names[item['id']] = next(x.split('=', 1)[1] for x in item['args'] if x.startswith('--job-name='))
for job in ids:
    if '--format=%A' in sys.argv:
        print(job.split('_', 1)[0])
    else:
        print(job + ' ' + names.get(job, 'unknown') + ' RUNNING 0:01 node')
''')
    scancel = directory / 'scancel'
    scancel.write_text('''#!/usr/bin/env python3
import os, pathlib, sys
if os.environ.get('FAIL_SCANCEL') == '1':
    print('deliberate scancel failure', file=sys.stderr)
    raise SystemExit(1)
root = pathlib.Path(os.environ['FAKE_SLURM_DIR'])
(root / 'cancelled').write_text('\\n'.join(sys.argv[1:]))
''')
    for p in (sbatch, scontrol, squeue, scancel):
        p.chmod(0o755)


@unittest.skipUnless(BWA.is_file(), 'bwa-mem2 is needed to build a reference index')
class SlurmWorkflow(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='rust-hic-slurm-test-')
        self.root = Path(self.temp.name)
        self.scratch = self.root / "scratch ' shared"
        self.scratch.mkdir()
        self.scheduler = self.root / 'scheduler'
        self.scheduler.mkdir()
        write_scheduler(self.scheduler)
        self.env = dict(os.environ, PATH=str(self.scheduler) + os.pathsep + os.environ['PATH'],
                        FAKE_SLURM_DIR=str(self.scheduler))
        self.reference = self.root / 'reference.fa'
        shutil.copyfile(DATA / 'genome.fa', self.reference)
        result = call(BWA, 'index', self.reference)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.base = {
            'assembly': 'test', 'chromsizes': str(DATA / 'chrom.sizes'),
            'scratch_dir': str(self.scratch), 'keep_pairs': True,
            'resolutions': [1000, 2000], 'balance_max_iterations': 2,
            'resources': {'sort_jobs': 2, 'sort_threads': 1, 'sort_run_mib': 1,
                          'merge_workers': 2, 'dedup_workers': 2, 'matrix_threads': 2},
            'fastq': {'read1': str(DATA / 'R1.fastq.gz'), 'read2': str(DATA / 'R2.fastq.gz'),
                      'reference': str(self.reference), 'aligner_index': str(self.reference),
                      'chunks': 3, 'jobs': 2, 'threads_per_job': 2, 'index_threads_per_mate': 1},
        }

    def tearDown(self):
        self.temp.cleanup()

    def submit(self, output, env=None, extra=()):
        config = self.root / (output.name + '.json')
        config.write_text(json.dumps({**self.base, 'output_dir': str(output)}))
        return call(BIN, 'run', config, '--mode', 'slurm', '--account', 'test', '--partition', 'cpu',
                    '--time', '00:30:00', *extra, env=env or self.env)

    def stage(self, work, name, chunk=None):
        env = dict(self.env)
        if chunk is not None:
            env['SLURM_ARRAY_TASK_ID'] = str(chunk)
        return call('/bin/sh', work / f'slurm-{name}.sh', env=env)

    def test_run_modes_and_slurm_preflight(self):
        out = self.root / 'checked'
        config = self.root / 'checked.json'
        config.write_text(json.dumps({**self.base, 'output_dir': str(out)}))
        for mode in ((), ('--mode', 'local')):
            result = call(BIN, 'run', config, *mode, '--check', env=self.env)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('Preflight passed', result.stdout)
            self.assertFalse(out.exists())
        result = call(BIN, 'run', config, '--mode', 'slurm', '--check',
                      '--account', 'test', '--map-mem', '24G', env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Slurm preflight passed', result.stdout)
        self.assertIn('Chunks: 3 (up to 2 concurrent)', result.stdout)
        self.assertIn('map 24G', result.stdout)
        self.assertFalse(out.exists())
        self.assertFalse((self.scheduler / 'submissions.jsonl').exists())
        env = dict(self.env, FAIL_SCONTROL='version')
        result = call(BIN, 'run', config, '--mode', 'slurm', '--check', env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('scontrol --version failed', result.stderr)
        self.assertFalse(out.exists())
        result = call(BIN, 'run', config, '--account', 'test', env=self.env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Slurm options require --mode slurm', result.stderr)
        self.assertFalse(out.exists())

    def test_explicit_local_matches_default_local(self):
        outputs = []
        for name, mode in [('default-local', ()), ('explicit-local', ('--mode', 'local'))]:
            out = self.root / name
            config = self.root / f'{name}.json'
            config.write_text(json.dumps({**self.base, 'output_dir': str(out)}))
            result = call(BIN, 'run', config, *mode, env=self.env)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(out)
        self.assertEqual((outputs[0] / 'filtered-nodups.pairs').read_bytes(),
                         (outputs[1] / 'filtered-nodups.pairs').read_bytes())
        with h5py.File(outputs[0] / 'contacts.mcool') as a, h5py.File(outputs[1] / 'contacts.mcool') as b:
            for resolution in (1000, 2000):
                for col in ('bin1_id', 'bin2_id', 'count'):
                    np.testing.assert_array_equal(a[f'resolutions/{resolution}/pixels/{col}'][:],
                                                  b[f'resolutions/{resolution}/pixels/{col}'][:])

    def test_end_to_end_matches_local_and_retries_are_safe(self):
        out = self.root / "distributed result's %J"
        submitted = self.submit(out, extra=('--keep-bam', '--sbatch-arg=--qos=normal'))
        self.assertEqual(submitted.returncode, 0, submitted.stderr)
        status_line = next(line for line in submitted.stdout.splitlines() if line.strip().startswith('Status:'))
        self.assertEqual(shlex.split(status_line.strip().removeprefix('Status: ')),
                         ['rust-hic', 'run', 'status', str(out)])
        jobs = json.loads((out / 'slurm-jobs.json').read_text())
        group_id = jobs['group_id']
        self.assertEqual(len(group_id), 8)
        self.assertTrue(group_id.isascii() and group_id.isalnum())
        self.assertIn(f'Group ID: {group_id}', submitted.stdout)
        self.assertEqual((jobs['index'], jobs['alignment'], jobs['finish']), ('1001', '1002', '1003'))
        requests = [json.loads(line) for line in (self.scheduler / 'submissions.jsonl').read_text().splitlines()]
        self.assertEqual(len(requests), 3)
        self.assertIn(f'--job-name=rust-hic-{group_id}-main', requests[0]['args'])
        self.assertIn(f'--job-name=rust-hic-{group_id}-mapping', requests[1]['args'])
        self.assertIn(f'--job-name=rust-hic-{group_id}-main', requests[2]['args'])
        self.assertIn('--hold', requests[1]['args'])
        self.assertNotIn('--hold', requests[0]['args'])
        self.assertNotIn('--hold', requests[2]['args'])
        controls = [json.loads(line) for line in (self.scheduler / 'control.jsonl').read_text().splitlines()]
        self.assertEqual(controls, [
            ['update', 'JobId=1002_0', f'JobName=rust-hic-{group_id}-mapping-0'],
            ['update', 'JobId=1002_1', f'JobName=rust-hic-{group_id}-mapping-1'],
            ['update', 'JobId=1002_2', f'JobName=rust-hic-{group_id}-mapping-2'],
            ['release', '1002'],
        ])
        self.assertIn('--cpus-per-task=3', requests[1]['args'])
        self.assertIn('--array=0-2%2', requests[1]['args'])
        self.assertIn('--dependency=afterok:1001', requests[1]['args'])
        self.assertIn('--dependency=afterany:1002', requests[2]['args'])
        self.assertIn('--no-requeue', requests[2]['args'])
        self.assertIn('--qos=normal', requests[0]['args'])
        self.assertTrue(any("result's %%J/slurm/" in arg for arg in requests[1]['args']))
        work = Path((out / 'work-dir.txt').read_text().strip())
        result = self.stage(work, 'index')
        self.assertEqual(result.returncode, 0, result.stderr)
        # Two attempts of the same chunk may run concurrently after a requeue.
        with ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(lambda i: self.stage(work, 'map', i), [0, 0, 1, 2]))
        for result in results:
            self.assertEqual(result.returncode, 0, result.stderr)
        result = self.stage(work, 'map', 1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('slurm_chunk_reused', result.stderr)
        (self.scheduler / 'active').write_text('1002_0\n1003\n')
        status = call(BIN, 'run', 'status', out, env=self.env)
        self.assertEqual(status.returncode, 0, status.stderr)
        self.assertIn('Alignment chunks: 3/3', status.stdout)
        self.assertIn(f'Group ID: {group_id}', status.stdout)
        self.assertIn(f'1002_0 rust-hic-{group_id}-mapping-0 RUNNING', status.stdout)
        self.assertIn(f'1003 rust-hic-{group_id}-main RUNNING', status.stdout)
        finish = self.stage(work, 'finish')
        self.assertEqual(finish.returncode, 0, finish.stderr)
        self.assertFalse(work.exists())
        self.assertFalse((out / 'work-dir.txt').exists())
        self.assertTrue((out / 'alignments.bam').is_file())
        self.assertEqual(json.loads((out / 'status.json').read_text())['status'], 'complete')
        report = json.loads((out / 'report.json').read_text())
        self.assertEqual(report['input']['aligner'], 'bwa-mem3')
        self.assertEqual(report['counts']['input_read_pairs'], 1024)
        self.assertTrue(any('Slurm' in item['stage'] for item in report['timings']))
        self.assertIn('Chunk preparation (Slurm finish job)',
                      [item['stage'] for item in report['timings']])
        self.assertIn('Alignment through chunk preparation (includes waits)',
                      [item['stage'] for item in report['timings']])
        if SAMTOOLS.is_file():
            self.assertEqual(call(SAMTOOLS, 'quickcheck', out / 'alignments.bam').returncode, 0)
            self.assertGreaterEqual(int(call(SAMTOOLS, 'view', '-c', out / 'alignments.bam').stdout), 2048)
        local = self.root / 'local'
        config = self.root / 'local.json'
        config.write_text(json.dumps({**self.base, 'output_dir': str(local)}))
        result = call(BIN, 'run', config, '--keep-bam', env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((out / 'filtered-nodups.pairs').read_bytes(),
                         (local / 'filtered-nodups.pairs').read_bytes())
        with h5py.File(out / 'contacts.mcool') as a, h5py.File(local / 'contacts.mcool') as b:
            for resolution in (1000, 2000):
                for col in ('bin1_id', 'bin2_id', 'count'):
                    np.testing.assert_array_equal(a[f'resolutions/{resolution}/pixels/{col}'][:],
                                                  b[f'resolutions/{resolution}/pixels/{col}'][:])
                np.testing.assert_allclose(a[f'resolutions/{resolution}/bins/weight'][:],
                                           b[f'resolutions/{resolution}/bins/weight'][:], equal_nan=True)

    def test_missing_chunk_fails_without_publishing_matrix(self):
        out = self.root / 'missing'
        result = self.submit(out)
        self.assertEqual(result.returncode, 0, result.stderr)
        work = Path((out / 'work-dir.txt').read_text().strip())
        self.assertEqual(self.stage(work, 'index').returncode, 0)
        result = self.stage(work, 'map', 0)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.stage(work, 'finish')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('chunk 1 is incomplete', result.stderr)
        self.assertTrue((out / 'failure.json').is_file())
        self.assertFalse((out / 'contacts.mcool').exists())
        self.assertTrue(work.is_dir())
        status = call(BIN, 'run', 'status', out, env=self.env)
        self.assertIn('Workflow: failed', status.stdout)
        self.assertIn('Alignment chunks: 1/3', status.stdout)

    def test_index_retry_clears_only_its_failure(self):
        r2 = self.root / 'retry-r2.fastq.gz'
        shutil.copyfile(DATA / 'R2.fastq.gz', r2)
        self.base['fastq']['read2'] = str(r2)
        out = self.root / 'index-retry'
        self.assertEqual(self.submit(out).returncode, 0)
        work = Path((out / 'work-dir.txt').read_text().strip())

        unavailable = self.root / 'temporarily-unavailable.fastq.gz'
        r2.rename(unavailable)
        try:
            first = self.stage(work, 'index')
        finally:
            unavailable.rename(r2)
        self.assertNotEqual(first.returncode, 0)
        self.assertEqual(json.loads((out / 'failure.json').read_text())['phase'], 'index')

        retry = self.stage(work, 'index')
        self.assertEqual(retry.returncode, 0, retry.stderr)
        self.assertFalse((out / 'failure.json').exists())
        (self.scheduler / 'active').write_text('1002_0\n')
        status = call(BIN, 'run', 'status', out, env=self.env)
        self.assertIn('Workflow: active or waiting', status.stdout)
        self.assertNotIn('Workflow: failed', status.stdout)

        finish = self.stage(work, 'finish')
        self.assertNotEqual(finish.returncode, 0)
        self.assertEqual(json.loads((out / 'failure.json').read_text())['phase'], 'finish')
        self.assertEqual(self.stage(work, 'index').returncode, 0)
        self.assertEqual(json.loads((out / 'failure.json').read_text())['phase'], 'finish')

    def test_reused_index_and_restriction_filter_match_local(self):
        index = self.root / 'existing-index'
        result = call(BIN, 'index', DATA / 'R1.fastq.gz', DATA / 'R2.fastq.gz',
                      '--chunks', '3', '--threads-per-mate', '1', '-o', index, env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.base['fastq']['index_dir'] = str(index)
        del self.base['fastq']['chunks']  # Use the reusable index's chunk count.
        self.base['fragments'] = str(DATA / 'fragments.bed')
        self.base['compressed_pairs'] = True
        out = self.root / 'restriction'
        result = self.submit(out)
        self.assertEqual(result.returncode, 0, result.stderr)
        work = Path((out / 'work-dir.txt').read_text().strip())
        result = self.stage(work, 'index')
        self.assertEqual(result.returncode, 0, result.stderr)
        for chunk in range(3):
            result = self.stage(work, 'map', chunk)
            self.assertEqual(result.returncode, 0, result.stderr)
        result = self.stage(work, 'finish')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(index.is_dir())  # An externally owned index is retained.
        self.assertTrue((out / 'filtered-nodups.pairs.gz.px2').is_file())
        local = self.root / 'restriction-local'
        config = self.root / 'restriction-local.json'
        config.write_text(json.dumps({**self.base, 'output_dir': str(local)}))
        result = call(BIN, 'run', config, env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((out / 'filtered-nodups.pairs').read_bytes(),
                         (local / 'filtered-nodups.pairs').read_bytes())
        self.assertEqual(json.loads((out / 'filter-stats.json').read_text()),
                         json.loads((local / 'filter-stats.json').read_text()))
        with h5py.File(out / 'contacts.mcool') as a, h5py.File(local / 'contacts.mcool') as b:
            for col in ('bin1_id', 'bin2_id', 'count'):
                np.testing.assert_array_equal(a[f'resolutions/1000/pixels/{col}'][:],
                                              b[f'resolutions/1000/pixels/{col}'][:])

    def test_chunk_count_estimate_and_expected_pair_count(self):
        del self.base['fastq']['chunks']
        largest = max((DATA / f'R{i}.fastq.gz').stat().st_size for i in (1, 2))
        self.base['fastq']['chunk_size_bytes'] = (largest + 2) // 3
        self.base['fastq']['expected_pairs'] = 1025  # The data actually contains 1024.
        out = self.root / 'count'
        result = self.submit(out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((out / 'slurm-jobs.json').read_text())['chunks'], 3)
        work = Path((out / 'work-dir.txt').read_text().strip())
        self.assertEqual(self.stage(work, 'index').returncode, 0)
        for chunk in range(3):
            result = self.stage(work, 'map', chunk)
            self.assertEqual(result.returncode, 0, result.stderr)
        result = self.stage(work, 'finish')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('read 1024 pairs but expected 1025', result.stderr)
        self.assertFalse((out / 'contacts.mcool').exists())

    def test_changed_chunk_is_not_reused(self):
        out = self.root / 'changed-chunk'
        self.assertEqual(self.submit(out).returncode, 0)
        work = Path((out / 'work-dir.txt').read_text().strip())
        self.assertEqual(self.stage(work, 'index').returncode, 0)
        self.assertEqual(self.stage(work, 'map', 0).returncode, 0)
        pair_file = work / 'chunks/0000/align-0000/alignments.pairs'
        with pair_file.open('ab') as f:
            f.write(b'changed\n')
        result = self.stage(work, 'map', 0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('pairs file is missing or changed', result.stderr)

    @unittest.skipUnless(CHROMAP.is_file(), 'chromap is needed for this test')
    def test_chromap_chunks_match_local_and_preserve_bam(self):
        index = self.root / 'chromap.index'
        result = call(CHROMAP, '-i', '-r', self.reference, '-o', index)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.base['fastq'].update({'aligner': 'chromap', 'aligner_executable': str(CHROMAP),
                                   'aligner_index': str(index)})
        out = self.root / 'chromap'
        result = self.submit(out, extra=('--keep-bam',))
        self.assertEqual(result.returncode, 0, result.stderr)
        work = Path((out / 'work-dir.txt').read_text().strip())
        result = self.stage(work, 'index')
        self.assertEqual(result.returncode, 0, result.stderr)
        with ThreadPoolExecutor(max_workers=3) as pool:
            results = list(pool.map(lambda i: self.stage(work, 'map', i), range(3)))
        for result in results:
            self.assertEqual(result.returncode, 0, result.stderr)
        result = self.stage(work, 'finish')
        self.assertEqual(result.returncode, 0, result.stderr)
        if SAMTOOLS.is_file():
            self.assertEqual(call(SAMTOOLS, 'quickcheck', out / 'alignments.bam').returncode, 0)
        local = self.root / 'chromap-local'
        config = self.root / 'chromap-local.json'
        config.write_text(json.dumps({**self.base, 'output_dir': str(local)}))
        result = call(BIN, 'run', config, env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((out / 'filtered-nodups.pairs').read_bytes(),
                         (local / 'filtered-nodups.pairs').read_bytes())
        report = json.loads((out / 'report.json').read_text())
        self.assertEqual(report['input']['aligner'], 'chromap')
        self.assertEqual(report['counts']['input_read_pairs'], 1024)

    def test_modified_config_or_index_is_rejected(self):
        out = self.root / 'modified'
        self.assertEqual(self.submit(out).returncode, 0)
        work = Path((out / 'work-dir.txt').read_text().strip())
        self.assertEqual(self.stage(work, 'index').returncode, 0)
        manifest = work / 'fastq-index' / 'manifest.json'
        manifest.write_text(manifest.read_text() + '\n')
        result = self.stage(work, 'map', 0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('manifest changed', result.stderr)
        config = out / 'config.resolved.yaml'
        config.write_text(config.read_text() + '\n')
        result = self.stage(work, 'map', 0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('configuration changed', result.stderr)
        self.assertFalse((out / 'contacts.mcool').exists())

    def test_submission_failure_cancels_own_jobs_and_cancel_only_targets_active(self):
        out = self.root / 'submit-failure'
        env = dict(self.env, FAIL_SBATCH='slurm-map.sh')
        result = self.submit(out, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('deliberate scheduler failure', result.stderr)
        self.assertEqual((self.scheduler / 'cancelled').read_text(), '1001')
        self.assertTrue((out / 'failure.json').exists())
        other = self.root / 'cancel'
        result = self.submit(other)
        self.assertEqual(result.returncode, 0, result.stderr)
        jobs = json.loads((other / 'slurm-jobs.json').read_text())
        (self.scheduler / 'active').write_text(jobs['alignment'] + '\n' + jobs['finish'] + '\n999999\n')
        result = call(BIN, 'run', 'cancel', other, env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(set((self.scheduler / 'cancelled').read_text().splitlines()),
                         {jobs['alignment'], jobs['finish']})
        self.assertTrue(other.is_dir())

    def test_separate_runs_use_distinct_group_ids_and_old_status_still_works(self):
        groups = []
        for name in ('group-one', 'group-two'):
            out = self.root / name
            result = self.submit(out)
            self.assertEqual(result.returncode, 0, result.stderr)
            jobs_path = out / 'slurm-jobs.json'
            jobs = json.loads(jobs_path.read_text())
            groups.append(jobs['group_id'])
        self.assertNotEqual(groups[0], groups[1])
        # Older job records lack the new group ID and remain inspectable.
        jobs.pop('group_id')
        jobs_path.write_text(json.dumps(jobs))
        result = call(BIN, 'run', 'status', out, env=self.env)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_managed_scheduler_options_cannot_be_overridden(self):
        for option in ('--sbatch-arg=--dependency=afterok:99', '--sbatch-arg=--hold=yes'):
            with self.subTest(option=option):
                out = self.root / 'bad-extra'
                result = self.submit(out, extra=(option,))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('managed Slurm option', result.stderr)
                self.assertFalse(out.exists())

    def test_name_or_release_failure_cancels_all_submitted_jobs(self):
        for operation in ('update', 'release'):
            with self.subTest(operation=operation):
                out = self.root / f'control-{operation}'
                result = self.submit(out, env=dict(self.env, FAIL_SCONTROL=operation))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('deliberate scontrol failure', result.stderr)
                jobs = json.loads((out / 'slurm-jobs.json').read_text())
                expected = {job for job in (jobs['index'], jobs['alignment'], jobs['finish']) if job}
                self.assertEqual(set((self.scheduler / 'cancelled').read_text().splitlines()), expected)
                self.assertEqual(len(expected), 2 if operation == 'update' else 3)
                self.assertTrue((out / 'failure.json').is_file())

    def test_failed_rollback_reports_job_ids(self):
        out = self.root / 'cancel-failure'
        env = dict(self.env, FAIL_SCONTROL='update', FAIL_SCANCEL='1')
        result = self.submit(out, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('unable to cancel submitted Slurm jobs 1001,1002', result.stderr)
        self.assertIn('deliberate scancel failure', result.stderr)
        self.assertIn('deliberate scontrol failure', result.stderr)
        failure = json.loads((out / 'failure.json').read_text())
        self.assertIn('unable to cancel submitted Slurm jobs 1001,1002', failure['error'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
