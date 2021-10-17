#!/usr/bin/env python3
from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def fasta_to_dict(path):
    return {record.id: record for record in SeqIO.parse(path, 'fasta')}

print('Loading FASTA...')
ref = fasta_to_dict(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa')
variants = {'ref': ref}
for x in (Path(__file__).parent / 'aligned_to_reference').glob('*.fa'):
    if 'ref' in x.name:
        continue
    seqs = fasta_to_dict(x)
    print(f'{x} -> {list(seqs.keys())}')
    variants[x.name[:-4]] = seqs

# For each chromossome, run multiple sequence alignement
for chromossome in ref:
    print(f'Multiple Sequence Alignment on {chromossome}')
    records = [SeqRecord(seqs[chromossome].seq, id=variant) for variant, seqs in variants.items()]
    chromossome_file = Path(__file__).parent / f'aligned_to_reference/{chromossome}.fa'
    SeqIO.write(records, chromossome_file, 'fasta')
    subprocess.run(['muscle', '-in', chromossome_file, '-out', f'{chromossome_file}_aligned.fasta', '-maxiters', '1', '-diags1', '-sv'], check=True)
