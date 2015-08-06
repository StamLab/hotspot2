#!/bin/bash

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

set -x -e -o pipefail

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

usage="Usage: simulate.tags.sh [--help] [--contig=] [--tmpdir=] [--seed=] <fragments> <uniq-mapability> <output file>"

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
echo "PARAM:uniquely_mapping_file:$uniq_mapping_file"
echo "PARAM:outfile:$outfile"
echo "PARAM:tmpdir:$tmpdir"
echo "PARAM:seed:$seed"

######
# Theorectically speaking, nothing needs to be
# changed below.
######

echo "BEGIN:Compute distributions. (`date -u`)"

FRAGCOUNTS=${tmpdir}/fragcount
FRAGLENS=${tmpdir}/fraglens
UNIQCOUNTS=${tmpdir}/uniqcounts

bedops $read_opts -u $frags \
	| awk \
		-v OFS="\t" \
		-v fraglens=${FRAGLENS} \
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
> ${FRAGCOUNTS}

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
> ${UNIQCOUNTS}

echo "END:Compute distributions. (`date -u`)"

files=""

# iterate contigs

while read chr frag uniq; do

echo "chr: $chr frag: $frag uniq: $uniq" 

    # sample line numbers and lengths
    LINENUMS=${tmpdir}/linenums.$chr
    TAGCOORDS=${tmpdir}/tags.${chr}.bed
    TAGCOUNTS=${tmpdir}/tagcounts.${chr}.bed

echo creating line numbers ${LINENUMS} 
    python $SCRIPT_DIR/subscripts/sample_linenums.py ${FRAGLENS} $uniq $frag $seed \
      | sort -k1,1g --buffer-size=4G \
      > ${LINENUMS}

if [ ! -s ${LINENUMS} ]; then
  echo Could not create ${LINENUMS}
  exit 1
else
  echo finished ${LINENUMS} 
fi
	# convert linenums and lengths to coordinates
echo creating tag coordinates ${TAGCOORDS}
	bedops --chrom $chr -u ${uniq_mapping_file} \
		| python $SCRIPT_DIR/subscripts/convert_linenums_to_coords.py ${LINENUMS} \
		| sort-bed - \
	> ${TAGCOORDS}

if [ ! -s ${TAGCOORDS} ]; then
  echo Could not create ${TAGCOORDS}
  exit 1
else
  echo finished ${TAGCOORDS}
fi

	# make a tag counts file
echo creating tag counts ${TAGCOUNTS}
	bedops -m ${TAGCOORDS} | bedops --chop 1 - \
		| bedmap --faster --delim "\t" --echo --count - ${TAGCOORDS} \
		| awk -v OFS="\t" '{ print $1, $2, $3, "i", $4; }' \
	> ${TAGCOUNTS}

if [ ! -s ${TAGCOUNTS} ]; then
  echo Could not create ${TAGCOUNTS}
  exit 1
else
  echo finished ${TAGCOUNTS}
fi

	files="$files ${TAGCOUNTS}"

#end while loop

echo FINISHED $chr
echo "files: $files"
done < <(join -j 1 ${FRAGCOUNTS} ${UNIQCOUNTS})

# merge together

bedops -u $files > ${tmpdir}/tagcounts.bed

# if starch, else copy

outfilename=$(basename "$outfile")
extension="${outfilename##*.}"

if [ "$extension" == "starch" ]; then
	starch ${tmpdir}/tagcounts.bed > $outfile
else
	rsync ${tmpdir}/tagcounts.bed $outfile
fi
