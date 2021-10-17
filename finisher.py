#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
import sam
import subprocess

from Bio import SeqIO

genome_paths = list((Path(__file__).parent / 'draft').glob('*_yeast_genoma.fa'))
country_genome = {x.name[:-16]: x for x in genome_paths}
print(f"{len(country_genome)} genomes found: {', '.join(country_genome)}")


# Align to reference genome with minimap2
def minimap2_to_sam(reference, draft, output):
    with open(output, 'w') as f:
        subprocess.run(['minimap2', '-cax', 'asm5', '-t24', reference, draft], check=True, stdout=f)


def minimap2_to_paf(reference, draft, output):
    with open(output, 'w') as f:
        subprocess.run(['minimap2', '-cx', 'asm5', '--cs', '-t24', reference, draft], check=True, stdout=f)


def ensure_sam():
    # Ensure SAM files were created by minimap2
    for country, path in country_genome.items():
        sam = Path(__file__).parent / f'aligned_to_reference/{country}.sam'
        paf = Path(__file__).parent / f'aligned_to_reference/{country}.paf'
        if sam.exists():
            continue
        seqs = list(SeqIO.parse(path, 'fasta'))
        print(f'{country} - {len(seqs)} entries')
        sizes = set()
        for seq in seqs:
            size = len(seq)
            sizes.add(size)
        print(f'Entry sizes: {sorted(sizes)}')
        print(f'Mapping {country} with minimap2...')
        minimap2_to_sam('S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', path, sam)
        minimap2_to_paf('S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', path, paf)


# def counts(s):
#     return {c: s.count(c) for c in set(s)}
# 
# print(largest, [str(largest[1]).count(c) for c in ('A', 'T', 'G', 'C', 'N')])


def consensus(sam_path):
    mapped_len = 0
    count = 0
    chromosome_entries = defaultdict(list)
    for entry in sam.SamReader(sam_path):
        reference_len = sam.cigar_reference_length(entry.cigar)
        chromosome_entries[entry.rname].append((entry.pos, reference_len, entry))
        mapped_len += reference_len
        count += 1
    print(f'Total length of mapped segments: {mapped_len} (average length: {mapped_len / count})')
    mapped_len = 0
    for chromosome, entries in chromosome_entries.items():
        print(chromosome)
        entries.sort()
        # print(entries)
        prev_end = 0
        for entry in entries:
            end = entry[0] + entry[1]
            if prev_end > entry[0]:
                print(f'Overlap at {entry[0]}')
            mapped_len += end - max(prev_end, entry[0])
            prev_end = end
    print(f'Mapped length without overlaps: {mapped_len}')

if __name__ == '__main__':
    ensure_sam()
    for aligned in (Path(__file__).parent / 'aligned_to_reference').glob('*.sam'):
        consensus(aligned)
