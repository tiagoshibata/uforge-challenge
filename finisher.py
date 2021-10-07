from collections import defaultdict
import pprint
from pathlib import Path

from Bio import Align, SeqIO
import edlib

genome_paths = list((Path(__file__).parent / 'draft').glob('*_yeast_genoma.fa'))
country_genome = {x.name[:-16]: x for x in genome_paths}
print(f"{len(country_genome)} genomes found: {', '.join(country_genome)}")
largest = (0, None)

for country, path in country_genome.items():
    seqs = list(SeqIO.parse(path, 'fasta'))
    print(f'{country} - {len(seqs)} entries')
    # sizes = defaultdict(int)
    sizes = set()
    total_size = 0
    for seq in seqs:
        size = len(seq)
        sizes.add(size)
        total_size += size
        if size > largest[0]:
            largest = (size, seq)
    # print(sorted(sizes))
    print(total_size)
    # pprint.pprint(data)

print(largest)

ref = list(SeqIO.parse(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', 'fasta'))
for chromosome in ref:
    print(f'Aligning to {chromosome.description} ({len(chromosome.seq)})')
    # print(set(largest[1].seq))  # A, T, G, C, N
    # print(set(chromosome.seq))
    print(edlib.align(largest[1].seq[::-1], chromosome.seq, 'HW', 'locations', 100000))
    # aligner = Align.PairwiseAligner()
    # aligner.mode = 'local'
    # print(chromosome.seq, largest[1])
    # alignments = aligner.align(chromosome.seq, largest[1].seq)
    # print(alignments[0])
