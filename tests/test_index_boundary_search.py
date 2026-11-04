#!/usr/bin/env python3
"""Check FASTQ boundary search when the R2 estimate is wrong in either direction."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get('RUST_HIC_BIN', ROOT / 'target/release/rust-hic'))


def run(*args):
    return subprocess.run([str(arg) for arg in args], capture_output=True, timeout=120)


class IndexBoundarySearch(unittest.TestCase):
    def test_searches_both_sides_and_preserves_every_record(self):
        with tempfile.TemporaryDirectory(prefix='fastq-boundary-search-') as d:
            tmp = Path(d)
            read1 = tmp / 'r1.fastq.gz'
            count = 120_000
            first = bytearray()
            for i in range(count):
                first.extend(f'@read{i:06d}/1\n'.encode())
                first.extend(b'A' * 100 + b'\n+\n' + b'I' * 100 + b'\n')
            read1.write_bytes(gzip.compress(first, compresslevel=1))
            for first_length, sign in [(20, -1), (180, 1)]:
                with self.subTest(r2_first_length=first_length):
                    read2 = tmp / f'r2-{first_length}.fastq.gz'
                    mate = bytearray()
                    for i in range(count):
                        length = first_length if i < count // 2 else 200 - first_length
                        mate.extend(f'@read{i:06d}/2\n'.encode())
                        mate.extend(b'A' * length + b'\n+\n' + b'I' * length + b'\n')
                    read2.write_bytes(gzip.compress(mate, compresslevel=1))
                    index = tmp / f'index-{first_length}'
                    args = [BIN, 'index', read1, read2, '--chunks', '2',
                            '--threads-per-mate', '1']
                    result = run(*args, '--search-mib', '16', '-o', index)
                    self.assertEqual(result.returncode, 0, result.stderr.decode(errors='replace'))
                    manifest = json.loads((index / 'manifest.json').read_text())
                    boundary = manifest['boundaries'][1]
                    predicted = (manifest['inputs'][1]['uncompressed_size'] *
                                 boundary['r1_uncompressed_offset'] //
                                 manifest['inputs'][0]['uncompressed_size'])
                    delta = boundary['r2_uncompressed_offset'] - predicted
                    self.assertGreater(sign * delta, 4 * 1024 * 1024)
                    for number, expected in [(1, first), (2, mate)]:
                        chunks = []
                        for chunk in (0, 1):
                            result = run(BIN, 'index', '--internal', 'stream', index, str(chunk), str(number))
                            self.assertEqual(result.returncode, 0, result.stderr.decode(errors='replace'))
                            chunks.append(result.stdout)
                        self.assertEqual(b''.join(chunks), expected)
                    too_narrow = tmp / f'narrow-{first_length}'
                    result = run(*args, '--search-mib', '1', '-o', too_narrow)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(b'R2 ID not found', result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
