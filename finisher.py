import pprint

from Bio import SeqIO
from pathlib import Path

genome_paths = list((Path(__file__).parent / 'draft').glob('*_yeast_genoma.fa'))
country_genome = {x.name[:-16]: x for x in genome_paths}
print(f"{len(country_genome)} genomes found: {', '.join(country_genome)}")

for country, path in country_genome.items():
    data = list(SeqIO.parse(path, 'fasta'))
    print(f'{country} - {len(data)} entries')
    # pprint.pprint(data)
