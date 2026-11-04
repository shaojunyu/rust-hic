#!/usr/bin/env python3
"""Regenerate the small deterministic dataset used by tests/test_core.py."""
import gzip
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent
rng = random.Random(824)
chromosomes = {name: ''.join(rng.choices('ACGT', k=25000)) for name in ('chr1', 'chr2')}
(HERE / 'genome.fa').write_text(''.join(f'>{name}\n{seq}\n' for name, seq in chromosomes.items()))
(HERE / 'chrom.sizes').write_text(''.join(f'{name}\t{len(seq)}\n' for name, seq in chromosomes.items()))
(HERE / 'fragments.bed').write_text(''.join(f'{name}\t0\t10\n{name}\t10\t{len(seq)}\n' for name, seq in chromosomes.items()))
reads = {1: bytearray(), 2: bytearray()}
names = tuple(chromosomes)
for i in range(1024):
    for mate in (1, 2):
        chrom = names[(i + mate) % 2]
        pos = rng.randrange(20, 24800)
        seq = chromosomes[chrom][pos:pos + 100] if i < 1020 else 'N' * 100
        reads[mate].extend(f'@read{i:04}/{mate}\n{seq}\n+\n{"I" * 100}\n'.encode())
for mate in (1, 2):
    with (HERE / f'R{mate}.fastq.gz').open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as out:
            out.write(reads[mate])
header = '## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n'
a = 'same_fragment\tchr1\t1\tchr1\t2\t+\t-\tUU\n'
b = 'retained_cis\tchr1\t1\tchr1\t11\t+\t-\tUU\n'
c = 'retained_trans\tchr1\t20\tchr2\t30\t+\t+\tUU\n'
duplicate = 'duplicate_trans\tchr1\t20\tchr2\t30\t+\t+\tUU\n'
(HERE / 'chunk1.pairs').write_text(header+c+a)
(HERE / 'chunk2.pairs').write_text(header+duplicate+b)
print('Generated test data in', HERE)
