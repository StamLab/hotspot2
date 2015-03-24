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

usage="Usage: random.tags.sh [--help] [--tmpdir=] [--seed] <tags> <uniq-mapability> <output BED file>"

params=$(getopt -o '' -l tmpdir:,seed:,help -n "random.tags.sh" -- "$@")
eval set -- "$params"

while true; do
	case "$1" in
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
		--help) echo -e $usage; exit 0; ;;
		--) shift; break; ;;
		*) echo -e "\Fatal error!\n"; exit 1; ;;
	esac	
done

tags=$1
uniq_mapping_file=$2
outfile=$3

if [ $# -lt 3 ]; then
	echo "ERROR: Missing required arguments!"
	echo $usage
	exit 1
fi

if [ ! -r "$tags" ]; then
	echo "ERROR: Tag file cannot be read!"
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

echo "PARAM:tagfile:$tags"
echo "PARAM:unqiuely_mapping_file:$uniq_mapping_file"
echo "PARAM:outfile:$outfile"
echo "PARAM:tmpdir:$tmpdir"
echo "PARAM:seed:$seed"

######
# Theorectically speaking, nothing needs to be
# changed below.
######

set -o pipefail

######
# The scripts translates line numbers into
# genomic coordinates.
######

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


rm -f ${tmpdir}/cutcounts.bed

cat $tags \
	| awk -v OFS="\t" '{ c[$1] += $5; } \
		END { for(chr in c) { print chr, c[chr]; } }' \
	| sort -k1,1 \
> ${tmpdir}/tagcounts

while read chr n; do

	# Generate random line numbers and sort

	bedextract $chr $uniq_mapping_file \
		| awk -v OFS="\t" \
			-v n=$n \
			-v seed=$seed \
			'BEGIN { srand(seed); } \
			{ total += $3-$2; } \
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
	>> ${tmpdir}/cutcounts.bed

done < ${tmpdir}/tagcounts

rsync ${tmpdir}/cutcounts.bed $outfile

