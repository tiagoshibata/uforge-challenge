import random
import sys

with open(sys.argv[1]) as f:
    for l in f.readlines():
        if l[0] == '>':
            print(l, end='')
        else:
            l = l.rstrip()
            print('\n'.join(l[i:i + 80] for i in range(0, len(l), 80)))
