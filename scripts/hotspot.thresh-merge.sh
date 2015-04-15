#/bin/bash

set -o pipefail

######
# Hawtspaught (BEDOPS-enabled)
# Jeff Vierstra
# 3/2015
#
# Based on Shane Neph's original implementation
######

######
# The following parameters must be set for the 
# software to function properly.
######

usage="Usage: hotspot.thresh-merge.sh [--help] [--tmpdir=] [--q-thresh=0.01] <unthresholded hotspots> <output hotspots file>"

params=$(getopt -o '' -l q-value:,tmpdir:,help -n "hotspot.thresh-merge.sh" -- "$@")
eval set -- "$params"

while true; do
	case "$1" in
		--q-value) 
			case "$2" in
				"") echo "ERROR: No q-value specified!"; exit 1; ;;
				*) q_thresh=$2; shift 2; ;;
			esac ;;
		--tmpdir) 
			case "$2" in
				"") echo "ERROR: No TMPDIR specified!"; exit 1; ;;
				*) tmpdir=$2; shift 2; ;;
			esac ;;
		--help) echo -e $usage; exit 0; ;;
		--) shift; break; ;;
		*) echo "Fatal error!"; exit 1; ;;
	esac	
done

hotspots=$1
outfile=$2

if [ $# -lt 2 ]; then
	echo "ERROR: Missing required arguments!"
	echo $usage
	exit 1
fi

if [ ! -r "$hotspots" ]; then
	echo "ERROR: Hotspots file cannot be read!"
	exit 1
fi

# If no TMPDIR set make one
# Else test if we have permissions, etc. to make one

if [ -z "$tmpdir" ]; then
	tmpdir=$(mktemp -d)
else
	mkdir -p $tmpdir || { echo "ERROR: Cannot create TMPDIR!"; exit 1; }
fi

# If no q-value set default of 0.01

if [ -z "$q_thresh" ]; then
	q_thresh=0.01
fi

# q-value threshold

cat $hotspots | awk -v OFS="\t" -v q_thresh=$q_thresh \
		'$11 > q_thresh { next; } { print; }'\
> ${tmpdir}/thresholded.${q_thresh}.bed

# merge down

bedops --range 75 -u ${tmpdir}/thresholded.${q_thresh}.bed \
	| bedops -m -\
	| bedops --range -75 -u - \
> ${tmpdir}/thresholded.merged.${q_thresh}.bed
	
# return merged with highest z-score

cut -f1-4,10 ${tmpdir}/thresholded.${q_thresh}.bed \
	| bedmap --delim "\t" --echo --max ${tmpdir}/thresholded.merged.${q_thresh}.bed - \
	| awk -v OFS="\t" '{ print $1, $2, $3, "id-"NR, $4; }' \
> $outfile
