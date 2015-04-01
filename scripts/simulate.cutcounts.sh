#!/bin/bash

set -o pipefail

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

usage="Usage: simulate.tags.sh [--help] [--contig=] [--tmpdir=] [--seed=] [--starch-output] <fragments> <uniq-mapability> <output file>"

params=$(getopt -o '' -l contig:,tmpdir:,seed:,starch-output,help -n "simulate.tags.sh" -- "$@")
eval set -- "$params"

while true; do
	case "$1" in
		--contig) 
			case "$2" in
				"") echo "ERROR: No contig specified!"; exit 1; ;;
				*) read_opts="--chrom $2"; shift 2; ;;
			esac ;;
		--tmpdir) 
			case "$2" in
				"") echo "ERROR: No TMPDIR specified!"; exit 1; ;;
				*) tmpdir=$2; shift 2; ;;
			esac ;;
		--seed) 
			case "$2" in
				"") echo "ERROR: No seed specified!"; exit 1; ;;
				*) seed=$2; shift 2; ;;
			esac ;;
		--starch-output)
			starch_output=1; shift; ;;
		--help) echo -e $usage; exit 0; ;;
		--) shift; break; ;;
		*) echo "Fatal error!"; exit 1; ;;
	esac	
done

frags=$1
uniq_mapping_file=$2
outfile=$3

if [ $# -lt 3 ]; then
	echo "ERROR: Missing required arguments!"
	echo $usage
	exit 1
fi

if [ ! -r "$frags" ]; then
	echo "ERROR: Fragments file cannot be read!"
	exit 1
fi

if [ ! -r "$uniq_mapping_file" ]; then
	echo "ERROR: Mappability file cannot be read!"
	exit 1
fi

# If no TMPDIR set make one
# Else test if we have permissions, etc. to make one

if [ -z "$tmpdir" ]; then
	tmpdir=$(mktemp -d)
else
	mkdir -p $tmpdir || { echo "ERROR: Cannot create TMPDIR!"; exit 1; }
fi

if [ -z "$seed" ]; then
	seed=1
fi

echo "PARAM:fragfile:$frags"
echo "PARAM:unqiuely_mapping_file:$uniq_mapping_file"
echo "PARAM:outfile:$outfile"
echo "PARAM:tmpdir:$tmpdir"
echo "PARAM:seed:$seed"

######
# Theorectically speaking, nothing needs to be
# changed below.
######

echo "BEGIN:Compute distributions. (`date -u`)"

bedops $read_opts -u $frags \
	| awk \
		-v OFS="\t" \
		-v fraglens=${tmpdir}/fraglens \
		'BEGIN { t = 0; } \
		{ c[$1] += 1; s[$3-$2] += 1; } \
		END { \
			for(contig in c) { \
				print contig, c[contig]; \
			} \
			for(len in s) { \
				print len, s[len] > fraglens; \
			} \
		}' \
	| sort -k1,1 \
> ${tmpdir}/fragcounts

bedops $read_opts -u $uniq_mapping_file \
	| awk \
		-v OFS="\t" \
		'{ c[$1] += $3-$2; } \
		END { \
			for(contig in c) { \
				print contig, c[contig]; \
			} \
		}' \
	| sort -k1,1 \
> ${tmpdir}/uniqcounts

echo "END:Compute distributions. (`date -u`)"

files=""

# iterate contigs

while read c n t; do

	# sample line numbers and lengths

	python <(cat <<__SCRIPT__

import numpy as np

fraglens = open("${tmpdir}/fraglens")

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

np.random.seed(${seed})

sampled_linenums = iter(np.random.choice(${t}, size = ${n}))
sampled_fraglens = iter(np.random.choice(lens, size = ${n}, p = prbs))

while 1:

	try:

		linenum = next(sampled_linenums)
		fraglen = next(sampled_fraglens)

		print "%d\t%d" % (linenum + 1, fraglen)

	except:

		break

__SCRIPT__
		) \
		| sort -k1,1g --buffer-size=4G \
	> ${tmpdir}/linenums

	# convert linenums and lengths to coordinates

	bedops --chrom $c -u ${uniq_mapping_file} \
		| python <(cat <<__SCRIPT__

import sys

uniqs = sys.stdin
linenums = open("${tmpdir}/linenums")

pos = -1

for line in linenums:

	(i, l) = line.strip().split('\t')
	(i, l) = (int(i), int(l))

	# if the line number is greater than the current
	# uniquely mapping segment, get the next one.

	while i > pos:

		rline = uniqs.readline()
		if not rline:
			break

		(contig, start, end) = rline.strip().split('\t')
		(start, end) = (int(start), int(end))

		pos += end - start
	
	# sample a second tag from the fragment length distribution

	s = end - (pos - i + 1)

	print "%s\t%d\t%d" % (contig, s, s + 1)
	print "%s\t%d\t%d" % (contig, s + l, s + l + 1);

linenums.close()

__SCRIPT__
		) \
		| sort-bed --max-mem 4G - \
	> ${tmpdir}/tags.${c}.bed

	# make a tag counts file

	bedops -m ${tmpdir}/tags.${c}.bed | bedops --chop 1 - \
		| bedmap --faster --delim "\t" --echo --count - ${tmpdir}/tags.${c}.bed \
		| awk -v OFS="\t" '{ print $1, $2, $3, "i", $4; }' \
	> ${tmpdir}/tagcounts.${c}.bed

	files="$files ${tmpdir}/tagcounts.${c}.bed"

#end while loop

done < <(join -j 1 ${tmpdir}/fragcounts ${tmpdir}/uniqcounts)

# merge together

bedops -u $files > ${tmpdir}/tagcounts.bed

# if starch, else copy

if [ -n "$starch_output" ]; then
	starch ${tmpdir}/tagcounts.bed > $outfile
else
	rsync ${tmpdir}/tagcounts.bed $outfile
fi

