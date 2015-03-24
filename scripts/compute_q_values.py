#!/bin/env python

import sys

import numpy as np

z_null = []

for line in open(sys.argv[1]):

	fields = line.strip().split('\t')

	z_null.append( float(fields[9]) )

z_null = np.sort(z_null)

p = np.arange(1, len(z_null) + 1, dtype = float) / len(z_null)
p = p[::-1]

import bisect

for line in open(sys.argv[2]):

	fields = line.strip().split('\t')

	z = float(fields[9])

	i = bisect.bisect_left(z_null, z)

	q = 0 if i == len(z_null) else p[i]

	print "%s\t%0.4f" % (line.strip(), q)