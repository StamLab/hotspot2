#!/bin/bash

######
# Hawtspaught (BEDOPS-enabled)
# Jeff Vierstra
# 3/2015
#
# This script will calculate the
# signal portion of tags (SPOT)
# from a tagcounts file and and
# a hotspots file
######


usage="Usage: spot.sh [--help] [--contig=] <tag counts> <hotspots>"

params=$(getopt -o '' -l contig:,help -n "spot.sh" -- "$@")
eval set -- "$params"

while true; do
	case "$1" in
		--contig) 
			case "$2" in
				"") echo "ERROR: No contig specified!"; exit 1; ;;
				*) read_opts="--chrom $2"; shift 2; ;;
			esac ;;
		--help) echo -e $usage; exit 0; ;;
		--) shift; break; ;;
		*) echo "Fatal error!"; exit 1; ;;
	esac	
done

tags=$1
hotspots=$2

if [ $# -lt 2 ]; then
	echo "ERROR: Missing required arguments!"
	echo $usage
	exit 1
fi

if [ ! -r "$tags" ]; then
	echo "ERROR: Tag file cannot be read!"
	exit 1
fi

if [ ! -r "$hotspots" ]; then
	echo "ERROR: Hotspots file cannot be read!"
	exit 1
fi

n=`bedops $read_opts -e -1 $tags $hotspots | awk '{ t += $5; } END { print t; }'`
total=`bedops $read_opts -u $tags | awk '{ t += $5; } END { print t; }'`

echo "SPOT: $(echo "scale=4; $n/$total" | bc)"
