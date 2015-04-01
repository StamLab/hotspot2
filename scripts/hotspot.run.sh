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

# A cut counts files can be made with the
# following code:
#
#	bedops -m cuts.bed | bedops --chop 1 - 
#		| bedmap --faster --delim "\t" --echo --count - cuts.bed \
#		| awk -v OFS="\t" '{ print $1, $2, $3, "i", $4; }' \
#	> cutcounts.bed

# For proper FDR calibration a simulated tag file
# is needed. There is an accompanying utility:
# "random.tags.sh" that will make one for you.


usage="Usage: hotspot.run.sh [--help] [--tmpdir=] [--contig=] <tag counts> <uniq-mapability> <output hotspots file>"

params=$(getopt -o '' -l contig:,tmpdir:,help -n "hotspot.run.sh" -- "$@")
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
		--help) echo -e $usage; exit 0; ;;
		--) shift; break; ;;
		*) echo "Fatal error!"; exit 1; ;;
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

echo "PARAM:tagfile:$tags"
echo "PARAM:uniquely_mapping_file:$uniq_mapping_file"
echo "PARAM:outfile:$outfile"
echo "PARAM:tmpdir:$tmpdir"
echo "PARAM:contig:$contig"

######
# Theorectically speaking, nothing needs to be
# changed below.
######

######
# Hardcoded variables
######

npasses=2 					# Number of passes
z_thresh=2 					# Min. z score to be reported (*not FDR*)
merge_distance=75 			# Min. distance to merge hotspots

local_window_size_min=100	# Min. local window size
local_window_size_max=150	# Max. local window size
local_window_size_step=25	# Windows sizes are stepped by this amount
local_window_sd_thresh=3	# Number of SDs need to pass threshold

local_window_sizes=(`seq $local_window_size_min $local_window_size_step $local_window_size_max`)

background_window_size=25000	# Background window size

#TODO: implement badspot pipeline
touch ${tmpdir}/badspots.bed
badspots="${tmpdir}/badspots.bed" # need to talk to Bob aboout how to make this file

######
# 
######

excluded="" # regions to be removed that are accumulated with each "pass"
included="" # putative hotspots identified with each "pass"

# Named pipes that are reused throughout the script

background_window_counts=${tmpdir}/background_window.counts.pipe
local_window_counts=${tmpdir}/local_window.counts.pipe

# Set the read command by format and 
# whether a specific contig is desired

read_command="bedops $read_opts -u $tags"	

# Make the I/O faster by uncompressing
# the cutcounts file to a local temp

eval $read_command > ${tmpdir}/tagcounts.bed

######
# Main iteration loop
######

for iteration in $(seq 1 $npasses); do
	
	# Step 1: Look for local enrichment of tags in
	# small windows (200-300 bp)

	# Step 2: Merge small windows as necesary

	# Step 3: Re-calculate z-score for windows

	# Step 4: Threshold and merge the hotspots from
	# each pass (z > 2 and within 75 bp of each other)
	# for exclusion in the next pass

	raw_hotspots=${tmpdir}/raw.hotspots.pass-${iteration}.bed
	scored_hotspots=${tmpdir}/scored.hotspots.pass-${iteration}.bed
	merged_thresholded_hotspots=${tmpdir}/merged-thresholded.hotspots.pass-${iteration}.bed

	echo "BEGIN:pass-$iteration:Selecting and merging windows. (`date -u`)"

	# Get number of tags in background window
	# Get number of uniquely mapping positions in background window
	
	rm -f $background_window_counts; mkfifo $background_window_counts

	bedops -n -1 ${tmpdir}/tagcounts.bed $badspots $excluded \
		| bedmap --faster --delim "\t" --prec 0 --range $background_window_size --echo --sum - \
		| bedops --range $background_window_size -u - \
		| bedmap --faster --delim "\t" --echo --bases-uniq - $uniq_mapping_file \
		| cut -f6- \
	> $background_window_counts &

	local_window_counts_files=""

	# Iterate over possible window sizes

	for window_size in ${local_window_sizes[@]}; do

		rm -f ${tmpdir}/local_window_${window_size}.counts.pipe; mkfifo ${tmpdir}/local_window_${window_size}.counts.pipe
		
		bedops -n -1 ${tmpdir}/tagcounts.bed $badspots $excluded \
			| bedmap --faster --prec 0 --range ${window_size} --sum - \
		> ${tmpdir}/local_window_${window_size}.counts.pipe &

		local_window_counts_files="$local_window_counts_files ${tmpdir}/local_window_${window_size}.counts.pipe"

	done
	
	######
	# Step 1 and 2:
	######

	# At each tag position get:
	# --# positions that are unqiuely mappable in background window
	# --# tags in background window
	# --# tags in a the local windows
	# Then:
	# --Find positions that are at least 3 SDs above expected
	# --Select biggest window size (per position) that passes the above criteria
	# Finally merge windows:
	# --Center is selected by averaging windows, window size is also averaged
	# --Output BED: chr window_left window_right "i" pos

	bedops -n -1 ${tmpdir}/tagcounts.bed $badspots $excluded \
		| cut -f1-3 - \
		| paste - $background_window_counts $local_window_counts_files \
		| awk -v OFS="\t" \
			-v n_windows=${#local_window_sizes[@]} \
			-v window_min=$local_window_size_min \
			-v window_step=$local_window_size_step \
			-v window_sd_thresh=$local_window_sd_thresh \
			'{ \
				max_window_size = -1; \
				max_o = 0; \
				max_e = 0; \
				\
				for(i = 0; i < n_windows; i += 1) { \
					window_size = (i * window_step) + window_min; \
					full_window_size = window_size * 2; \
					\
					if(full_window_size >= $5 || $4 == 0) { continue; } \
					\
					o = $(i+6); \
					\
					p = full_window_size / $5; \
					e = p * $4; \
					e_sigma = sqrt(e * (1-p)); \
					\
					thresh = 1 + e + (e_sigma * window_sd_thresh); \
					\
					if(o > thresh) { \
						max_window_size = window_size; \
						max_o = o; max_e = e; \
					} \
				} \
				\
				if(max_window_size > 0) { \
					print $1, $2, $3, "i", max_window_size, max_e, max_o; \
				} \
			}' \
		| awk -v OFS="\t" \
			'function abs(x) { \
				return ((x < 0.0) ? -x : x) \
			} \
			\
			BEGIN { chr = ""; } \
			{ \
				if (chr == $1 && abs(start - $2) < ($5 * 2)) { \
					center += $2; w += $5; e += $6; o += $7; n += 1; \
				} else { \
					if(chr != "") { \
						print chr, int(center/n) - int(w/n), int(center/n) + int(w/n), "i", int(center/n), e/n, o/n, n; \
					} \
					chr = $1; start = $2; center = $2; w = $5; e = $6; o = $7; n = 1; \
				} \
			} \
			END { print chr, int(center/n) - int(w/n), int(center/n) + int(w/n), "i", int(center/n), e/n, o/n, n; }' \
	> $raw_hotspots

	# Wait for everything to finish up
	
	wait

	# Clean-up named pipes

	rm -f $background_window_counts
	rm -f $local_window_counts_files

	echo "END:pass-$iteration:Selecting and merging windows. (`date -u`)"

	######
	# Step 3: Go back and calculate z-score;
	######

	# Note that we just remove tags in badspots

	echo "BEGIN:pass-$iteration:Re-calculating Z-scores from all data. (`date -u`)"

	# Get number of tags in background window
	# Get number of uniquely mapping positions in background window

	rm -f $background_window_counts; mkfifo $background_window_counts

	bedops -n -1 ${tmpdir}/tagcounts.bed $badspots \
		| bedmap --faster --delim "\t" --prec 0 --range $background_window_size --echo --sum $raw_hotspots - \
		| bedops --range $background_window_size -u - \
		| bedmap --faster --delim "\t" --echo --bases-uniq - $uniq_mapping_file \
		| cut -f9- \
	> $background_window_counts &

	# Get # of tags in hotspots

	rm -f $local_window_counts; mkfifo $local_window_counts

	bedops -n -1 ${tmpdir}/tagcounts.bed $badspots \
		| bedmap --faster --delim "\t" --prec 0 --echo-map-range --sum $raw_hotspots - \
		| cut -f2- \
	> $local_window_counts &

	# Re-calculate z-score
	# Input BED: chr window_left window_right "i" pos uniq_pos bkg left right obs
	# Output BED: chr left right "i" pos window_size expect expect_sigma, observed, z

	cut -f1-5 $raw_hotspots \
		| paste - $background_window_counts $local_window_counts \
	 	| awk -v OFS="\t" \
	 	'NF == 10 { \
	 		window_size = $3-$2; \
	 		\
	 		if(window_size >= $7 || $6 == 0) { next; } \
	 		\
	 		o = $10;\
	 		\
	 		p = window_size / $7; \
	 		e = p * $6; \
	 		e_sigma = sqrt(e * (1-p)); \
	 		\
	 		z = (o - e) / e_sigma; \
	 		\
	 		print $1, $8, $9, "i", $5, window_size, e, e_sigma, o, z; \
	 	}' \
	 	| sort-bed - \
	 > $scored_hotspots

	 # Wait for everything to finish up
	
	wait

	# Clean-up named pipe
	
	rm -f $background_window_counts
	rm -f $local_window_counts

	echo "END:pass-$iteration:Re-calculating Z-scores from all data. (`date -u`)"

	######
	# Step 4: Make thresholded and merged exclude list
	######
	
	echo "BEGIN:pass-$iteration:Generating excluded regions. (`date -u`)"

	cat $scored_hotspots \
	 	| awk -v OFS="\t" -v min_width=$min_width -v z_thresh=$z_thresh \
			'$3-$2 < min_width { next; } $10 < z_thresh { next; } { print; }' \
		| bedops --range $merge_distance -u - | bedops -m - | bedops --range -$merge_distance -u - \
		| sort-bed - \
	> $merged_thresholded_hotspots

	echo "END:pass-$iteration:Generating excluded regions. (`date -u`)"

	excluded="$excluded $merged_thresholded_hotspots"
	included="$included $scored_hotspots"

done

echo "BEGIN:Combine and threshold passes. (`date -u`)"

bedops -u $included \
	| awk -v OFS="\t" -v min_width=$min_width -v z_thresh=$z_thresh \
			'$3-$2 < min_width { next; } $10 < z_thresh { next; } { print; }' \
> ${tmpdir}/combined.all-passes.hotspots.bed

echo "END:Combine and threshold passes. (`date -u`)"

rsync ${tmpdir}/combined.all-passes.hotspots.bed $outfile

