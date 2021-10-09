from collections import defaultdict, namedtuple
import re

MatchedSeq = namedtuple('MatchedSeq', ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual', 'optional'])


def optional(tag, type, value):
    if type in ('A', 'Z'):
        return tag, str(value)
    if type == 'f':
        return tag, float(value)
    if type == 'i':
        return tag, int(value)
    raise NotImplementedError()


def parse_optional(s):
    return dict(optional(*entry.split(':')) for entry in s.split())


def cigar_reference_length(cigar):
    len = 0
    for entry in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        if entry.group(2) in ('M', 'D', 'N', '=', 'X'):
            len += int(entry.group(1))
    assert(len)
    return len


def cigar_max_insertions(cigar):
    len = 0
    longest = 0
    for entry in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        if entry.group(2) in ('I', 'S'):
            len += int(entry.group(1))
            longest = max(longest, len)
        else:
            len = 0
    return longest


def fields_to_MatchedSeq(fields):
    # Convert fields from strings
    return MatchedSeq(fields[0], fields[1], fields[2], int(fields[3]), int(fields[4]), fields[5],
            fields[6], int(fields[7]), int(fields[8]), fields[9], fields[10], parse_optional(fields[11]))


class SamReader:
    def __init__(self, path):
        self.path = path

    def __iter__(self):
        prev = None
        is_duplicate = False
        with open(self.path) as f:
            for line in f.readlines():
                if line[0] == '@':
                    continue
                fields = line.split(maxsplit=11)
                if fields[2] == '*':
                    continue
                # Ignore sequences with multiple matches to the reference
                if prev:
                    if fields[0] != prev:
                        if not is_duplicate:
                            yield fields_to_MatchedSeq(prev)
                        is_duplicate = False
                    else:
                        is_duplicate = True
                prev = fields
                

def main(sam_path):
    mapped_len = 0
    count = 0
    chromosome_entries = defaultdict(list)
    for entry in SamReader(sam_path):
        reference_len = cigar_reference_length(entry.cigar)
        chromosome_entries[entry.rname].append((entry.pos, reference_len, entry.seq))
        mapped_len += reference_len
        count += 1
    print(f'Mapped length: {mapped_len} (average length: {mapped_len / count})')
    mapped_len = 0
    for chromosome, entries in chromosome_entries.items():
        print(chromosome)
        entries.sort()
        prev_end = 0
        print(entries)
        for i in range(len(entries)):
            entry = entries[i]
            end = entry[0] + entry[1]
            mapped_len += end - max(prev_end, entry[0])
            prev_end = end
            if i < len(entries) - 1 and entries[i][0] == entries[i + 1][0]:
                seq1 = entries[i][2]
                seq2 = entries[i + 1][2]
                print(seq1)
                print(seq2)
                print(sum(a != b for a, b in zip(seq1, seq2)), min(len(seq1), len(seq2)))
    print(f'Mapped length without overlaps: {mapped_len}')
    print(f'Maximum contiguous insertions to reference genome: {max(cigar_max_insertions((x.cigar)) for x in SamReader(sam_path))}')


if __name__ == '__main__':
    from pathlib import Path

    main('./aligned_to_reference/Russia.sam')
    # for p in Path('aligned_to_reference').glob('*.sam'):
    #     main(p)
