import numpy as np
import sys

fraglens = open(sys.argv[1])
uniq = int(sys.argv[2])
frag = int(sys.argv[3])
seed = int(sys.argv[4])

lens = []
prbs = []

for line in fraglens:

    (l, n) = line.strip().split('\t')
    (l, n) = (int(l), float(n))

    lens.append(l)
    prbs.append(n)

prbs = np.array(prbs)
prbs /= np.sum(prbs)

fraglens.close()

np.random.seed(seed)

sampled_linenums = iter(np.random.choice(uniq, size = frag))
sampled_fraglens = iter(np.random.choice(lens, size = frag, p = prbs))

while 1:

    try:

        linenum = next(sampled_linenums)
        fraglen = next(sampled_fraglens)

        print "%d\t%d" % (linenum + 1, fraglen)

    except:

        break
