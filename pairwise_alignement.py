#!/usr/bin/env python3
from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def fasta_to_dict(path):
    return {record.id: record for record in SeqIO.parse(path, 'fasta')}

print('Loading FASTA...')
ref = fasta_to_dict(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa')
# variants = {'ref': ref}
variants = {}
for x in (Path(__file__).parent / 'aligned_to_reference').glob('*.fa'):
    if 'ref' in x.name:
        continue
    seqs = fasta_to_dict(x)
    print(f'{x} -> {list(seqs.keys())}')
    variants[x.name[:-3]] = seqs

# For each variant pair, run pairwise alignement
for v1 in variants:
    for v2 in variants:
        if v2 <= v1:
            continue
        # for chromossome in ref:
        print(f'Aligning {v1} x {v2}')
        with open(f'diff/{v1}_x_{v2}.sam', 'w') as f:
            subprocess.run(['minimap2', '-a', Path(__file__).parent / f'aligned_to_reference/{v1}.fa', Path(__file__).parent / f'aligned_to_reference/{v2}.fa'],
                    check=True, stdout=f)
