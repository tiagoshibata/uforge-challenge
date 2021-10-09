from collections import defaultdict
from pathlib import Path
from multiprocessing.pool import ThreadPool
import sam
import subprocess

from Bio import Align, SeqIO

ref = list(SeqIO.parse(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', 'fasta'))

genome_paths = list((Path(__file__).parent / 'draft').glob('*_yeast_genoma.fa'))
country_genome = {x.name[:-16]: x for x in genome_paths}
print(f"{len(country_genome)} genomes found: {', '.join(country_genome)}")


# Align to reference genome with minimap2
def minimap2(reference, draft, output):
    with open(output, 'w') as f:
        subprocess.run(['minimap2', '-cax', 'asm5', '-t24', reference, draft], check=True, stdout=f)


def minimap2_paf(reference, draft, output):
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
        minimap2('S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', path, sam)
        minimap2_paf('S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', path, paf)


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
        # Seek consensus on overlapping regions
        # i = 0
        # consensus = []
        # while i < len(entries):
        #     overlapping = [entries[i]]
        #     pos, len_ref, entry = entries[i]
        #     end_ref = pos + len_ref
        #     while i + 1 < len(entries) and entries[i + 1][0] < end_ref:
        #         i += 1
        #         end_ref = max(end_ref, entries[i][0] + entries[i][1])
        #         overlapping.append(entries[i])
        #     if len(overlapping) > 1:
        #         print(f'Overlapping sequences: {len(overlapping)}')
        #     i += 1
    print(f'Mapped length without overlaps: {mapped_len}')

if __name__ == '__main__':
    ensure_sam()
    consensus('./aligned_to_reference/Brasil.sam')
