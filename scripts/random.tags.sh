#!/bin/bash
#/bin/bash

######
# Hawtspaught (BEDOPS-enabled)
# Jeff Vierstra
# 3/2015
#
# This is a utility to make a simulated
# cutcounts file.
######

######
# The following parameters must be set for the 
# software to function properly.
######

outfile=
tmpdir=/tmp
cutcounts=/home/jvierstra/proj/ss.dnase/data/K562-P5-20140207/align/cutcounts.bed
uniq_mapping_file=/data/vol7/annotations/data/hg19/hg19.K36.mappable_only.bed

######
# Theorectically speaking, nothing needs to be
# changed below.
######

set -o pipefail

cat <<SCRIPT > ${tmpdir}/translate.py

import sys

start = 0 # region start
end = 0 # region end

pos = 0 # current line

fn = open(sys.argv[1])

for line in fn:
	
	(i, n) = line.strip().split('\t')
	(i, n) = (int(i), int(n))

	while i > pos:

		rline = sys.stdin.readline()
		if not rline:
			break

		(chr, start, end) = rline.strip().split('\t')
		(start, end) = (int(start), int(end))

		pos += end - start

	print "%s\t%d\t%d\ti\t%d" % (chr, end - (pos - i + 1), end - (pos - i + 1) + 1, n);

fn.close()

SCRIPT

rm -f $outfile

cat $cutcounts \
	| awk -v OFS="\t" '{ c[$1] += $5; } \
		END { for(chr in c) { print chr, c[chr]; } }' \
	| sort -k1,1 \
> ${tmpdir}/tagcounts

while read chr n; do

	# Generate random line numbers and sort

	bedextract $chr $uniq_mapping_file \
		| awk -v OFS="\t" \
			-v n=$n \
			'{ total += $3-$2; } \
			END { \
				for(i = 0; i < n; i++) { \
					print int(rand() * total); \
				} \
			 }' \
		| sort -g \
		| uniq -c \
		| awk -v OFS="\t" '{ print $2, $1 }' \
	> ${tmpdir}/linenums

	# Convert the line numbers to genomic positions

	bedextract $chr $uniq_mapping_file \
		| python ${tmpdir}/translate.py ${tmpdir}/linenums \
	>> $outfile

done < ${tmpdir}/tagcounts
