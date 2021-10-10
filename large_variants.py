from collections import defaultdict
from pathlib import Path
from multiprocessing.pool import ThreadPool
import sam
import subprocess
import pprint

from Bio import Align, SeqIO

ref = list(SeqIO.parse(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', 'fasta'))

genome_paths = list((Path(__file__).parent / './aligned_to_reference').glob('*.vcf'))
country_genome = {x.name[:-4]: x for x in genome_paths}
print(f"{len(country_genome)} genomes found: {', '.join(country_genome)}")

if __name__ == '__main__':
    for country, genome in country_genome.items():
        print(country)
        large_variations = []
        with genome.open() as f:
            for l in f.readlines():
                if l[0] == '#':
                    continue
                fields = l.split()
                if len(fields[3]) < 5000 and len(fields[4]) < 5000:
                    continue
                large_variations.append((fields[0], fields[1], len(fields[3]), len(fields[4])))
        pprint.pprint(large_variations)
