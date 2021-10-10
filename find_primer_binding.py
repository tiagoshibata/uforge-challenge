from pathlib import Path
import pprint

from Bio import SeqIO
from Bio.Seq import Seq

ref = list(SeqIO.parse(Path(__file__).parent / 'S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa', 'fasta'))
variants = list((Path(__file__).parent / 'aligned_to_reference').glob('*.fa'))

SC1 = 'AACGGTGAGAGATTTCTGTGC'
SC2 = 'AGCTGGCAGTATTCCCACAG'

primers = [
    Seq(SC2),
    Seq(SC1).reverse_complement(),
]


def truncate_at_stop(s):
    for i in range(0, len(s) - 3, 3):
        if s[i:i + 3] in ('TAG', 'TAA', 'TGA'):
            return s[:i + 3]


def binding_site(seq):
    sc1_bases = sc2_bases = None
    print(seq)
    f = list(SeqIO.parse(seq, 'fasta'))
    for primer in primers:
        assert(sum(primer in chromosome for chromosome in f) == 1)
        needle = str(primer)
        for chromosome in f:
            if primer not in chromosome:
                continue
            haystack = str(chromosome.seq)
            if needle == SC2:
                start = haystack.find(needle) + len(needle)
                sc2_bases = haystack[start:]
            else:
                start = haystack.find(needle)
                sc1_bases = haystack[:start:-1]

    assert(sc1_bases)
    assert(sc2_bases)
    
    return truncate_at_stop(sc1_bases), truncate_at_stop(sc2_bases)


if __name__ == '__main__':
    bases = {var: binding_site(var) for var in variants}
    pprint.pprint(bases)
    values = list(bases.values())
    print(all(x == values[0] for x in values))
