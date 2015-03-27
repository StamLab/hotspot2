#!/bin/bash

######
# Hawtspaught (BEDOPS-enabled)
# Jeff Vierstra
# 3/2015
#
# This script will downsample a cutcounts file 
# to the amount desired.
######

usage="Usage: downsample.sh [--help] [--tmpdir=] [--seed] <N> <tags> <output BED file>"

params=$(getopt -o '' -l tmpdir:,seed:,help -n "downsample.sh" -- "$@")
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
		*) echo "Fatal error!"; exit 1; ;;
	esac	
done

n=$1
tags=$2
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

######
# Theorectically speaking, nothing needs to be
# changed below.
######

set -o pipefail

# Count the total tags in the file

ntags=`bedops -u $tags | awk '{ n += $5; } END { print n; }'`

if [ ! "$n" -le "$ntags" ]; then
	echo "ERROR: Number of tags desired ($n) is greater or equal to total tags ($ntags)"
	exit 1
fi

# Sequential sample lines from the "uncompressed" file.

bedops -u $tags \
	| awk '{ for(i = 0; i < $5; i++) { print $1, $2, $3; } }' \
	| random-lines --num=$n --max=$ntags --seed=$seed \
> ${tmpdir}/cuts

# Package it back up into a per-position count file

bedops -m ${tmpdir}/cuts | bedops --chop 1 - \
	| bedmap --faster --delim "\t" --echo --count - ${tmpdir}/cuts \
	| awk -v OFS="\t" '{ print $1, $2, $3, "i", $4; }' \
> $outfile
