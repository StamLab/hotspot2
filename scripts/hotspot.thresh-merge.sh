q_thresh=$1
hotspots=$2
outfile=$3

tmpdir=/tmp

cat $hotspots | awk -v OFS="\t" -v q_thresh=$q_thresh \
		'$11 > q_thresh { next; } { print; }'\
> ${tmpdir}/thresholded.${q_thresh}.bed

bedops --range 75 -u ${tmpdir}/thresholded.${q_thresh}.bed \
	| bedops -m -\
	| bedops --range -75 -u - \
> ${tmpdir}/thresholded.merged.${q_thresh}.bed
	# return merged with highest z-score

cut -f1-4,10 ${tmpdir}/thresholded.${q_thresh}.bed \
	| bedmap --delim "\t" --echo --max ${tmpdir}/thresholded.merged.${q_thresh}.bed - \
> $outfile
