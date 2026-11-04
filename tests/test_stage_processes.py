#!/usr/bin/env python3
"""Exercise nested sort/merge stage processes with multiple runs and ranges."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))
HEADER = '## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n'


def run(*args):
    result = subprocess.run([str(a) for a in args], capture_output=True, text=True, timeout=120)
    if result.returncode:
        raise AssertionError(f'{args}: {result.stderr}')
    return result


class StageProcesses(unittest.TestCase):
    def test_sort_runs_and_parallel_merge(self):
        with tempfile.TemporaryDirectory(prefix='stage-processes-') as work:
            work = Path(work)
            sort = work / 'sort_pairs'
            merge = work / 'merge_pairs'
            sort.symlink_to(BIN.resolve())
            merge.symlink_to(BIN.resolve())
            chromsizes = work / 'chrom.sizes'
            chromsizes.write_text('chr1\t1000000\nchr2\t1000000\n')
            sorted_paths = []
            for source in range(2):
                path = work / f'input{source}.pairs'
                with path.open('w') as f:
                    f.write(HEADER)
                    for i in range(39999, -1, -1):
                        f.write(f'read{source}_{i:06}\tchr1\t{i * 2 + source + 1}\tchr2\t10\t+\t-\tUU\n')
                output = work / f'sorted{source}.pairs'
                result = run(sort, chromsizes, 'test', output, work, '2', '1', path)
                self.assertIn('runs=2', result.stderr)
                sorted_paths.append(output)
            result = run(merge, '4', work, *sorted_paths)
            self.assertIn('segments=4', result.stderr)
            rows = [line.split('\t') for line in result.stdout.splitlines() if not line.startswith('#')]
            self.assertEqual(len(rows), 80000)
            self.assertEqual([int(r[2]) for r in rows], list(range(1, 80001)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
