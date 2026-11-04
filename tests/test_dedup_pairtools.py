#!/usr/bin/env python3
"""Differential checks against pairtools' default scipy dedup backend."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))
PAIRTOOLS = os.environ.get('RUST_HIC_PAIRTOOLS', 'pairtools')
HEADER = '## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n'


def row(name, p1, p2, s1='+', s2='-'):
    return f'{name}\tchr1\t{p1}\tchr1\t{p2}\t{s1}\t{s2}\tUU\n'


def run(args, data):
    result = subprocess.run(args, input=data, text=True, capture_output=True, timeout=120)
    if result.returncode:
        raise AssertionError(f'{args}: {result.stderr}')
    return [x for x in result.stdout.splitlines() if not x.startswith('#')]


class PairtoolsDedup(unittest.TestCase):
    def test_fragment_layout_and_malformed_rows(self):
        fragment_header = HEADER.replace(
            'strand1 strand2 pair_type\n',
            'strand1 strand2 pair_type rfrag1 rfrag_start1 rfrag_end1 '
            'rfrag2 rfrag_start2 rfrag_end2\n',
        )
        def annotated(name, p1, p2):
            return row(name, p1, p2).rstrip('\n') + '\t1\t0\t100\t2\t100\t200\n'

        data = fragment_header + annotated('a', 10, 100) + annotated('b', 12, 102)
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / 'dedup_pairs'
            stage.symlink_to(BIN.resolve())
            expected = run([PAIRTOOLS, 'dedup'], data)
            for threads in (1, 4):
                self.assertEqual(run([str(stage), str(threads), '3'], data), expected)
                bad = subprocess.run([str(stage), str(threads), '3'],
                                     input=fragment_header + row('bad', 10, 100),
                                     text=True, capture_output=True, timeout=30)
                self.assertNotEqual(bad.returncode, 0)
                self.assertIn('pairs row has 8 fields', bad.stderr)

    def test_clean_default_and_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / 'input.pairs'
            source.write_text(HEADER + row('a', 10, 100) + row('b', 12, 102))
            sizes = tmp / 'chrom.sizes'
            sizes.write_text('chr1\t10000\n')
            for name, extra, count in [('default', [], 1), ('one', ['--max-mismatch', '1'], 2)]:
                output = tmp / f'{name}.pairs'
                args = [str(BIN), 'clean', str(source), '--chromsizes', str(sizes),
                        '--threads', '1', *extra, '-o', str(output)]
                result = subprocess.run(args, capture_output=True, text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(len([x for x in output.read_text().splitlines()
                                      if not x.startswith('#')]), count)
            config = tmp / 'run.yaml'
            config.write_text(f'assembly: test\nchromsizes: {sizes}\ninput_pairs:\n'
                              f'  - {source}\noutput_dir: {tmp / "run-result"}\n'
                              f'scratch_dir: {tmp}\ndedup_max_mismatch: 1\n'
                              'resolutions: [10]\nbalance_max_iterations: 1\n')
            result = subprocess.run([str(BIN), 'run', str(config)], capture_output=True,
                                    text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stderr)
            result_pairs = (tmp / 'run-result' / 'filtered-nodups.pairs').read_text()
            self.assertEqual(len([x for x in result_pairs.splitlines()
                                  if not x.startswith('#')]), 2)

    def test_mismatch_transitivity_and_chunk_boundary(self):
        rows = [row('a', 10, 100), row('b', 13, 103), row('c', 16, 106),
                row('outside', 20, 106), row('strand', 20, 106, '-', '-')]
        # Put a transitive duplicate chain across pairtools' 10,000-row boundary.
        rows += [row(f'filler{i}', 1000 + 10 * i, 5000) for i in range(9994)]
        rows += [row('edge-a', 1000000, 2000000), row('edge-b', 1000003, 2000003),
                 row('edge-c', 1000006, 2000006), row('edge-other', 1000020, 2000000)]
        data = HEADER + ''.join(rows)
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp) / 'dedup_pairs'
            stage.symlink_to(BIN.resolve())
            for mismatch in (0, 1, 3, 4, 17):
                with self.subTest(mismatch=mismatch):
                    actual = run([str(stage), '2', str(mismatch)], data)
                    expected = run([PAIRTOOLS, 'dedup', '--max-mismatch', str(mismatch)], data)
                    self.assertEqual(actual, expected)
            self.assertEqual(run([str(stage), '2', '3'], data), run([PAIRTOOLS, 'dedup'], data))


if __name__ == '__main__':
    unittest.main(verbosity=2)
