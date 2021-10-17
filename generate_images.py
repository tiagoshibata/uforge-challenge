#!/usr/bin/env python3
from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import sam

print('Loading FASTA...')
ref = Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa'
chromossomes = [r.id for r in SeqIO.parse(ref, 'fasta')]

for comp in (Path(__file__).parent / 'diff').glob('*.sam'):
    print(comp)
    data = sam.SamReader(comp)

