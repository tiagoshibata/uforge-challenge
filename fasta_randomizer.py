import random
import string
import sys

with open(sys.argv[1]) as f:
    for l in f.readlines():
        if l[0] == '>':
            print('>' + ''.join(random.choices(string.ascii_letters + string.digits, k=24)))
        else:
            l = l.rstrip()
            print('\n'.join(l[i:i + 80] for i in range(0, len(l), 80)))
