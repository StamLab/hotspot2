#!/bin/tcsh -ef
# author : Shane Neph
# date   : 2016
# proj   : Reimplement hotspot using BEDOPS; see bottom of file for improvements over original hotspot
#            - output includes badspot calculations and follow-on thresholded (by z-score) hotspot calculations
#            - may be called on real or simulated data.


################
# some defaults
################
set db = hg19
set datadir = ../data/$db
set basetmpdir = /tmp/`whoami`/hotspot
@ defaulttmpdir = 1
if ( $?TMPDIR ) then
  # use the environmental variable
  set basetmpdir = $TMPDIR
  @ defaulttmpdir = 0
endif
set tmpdir = $basetmpdir/$$


#########
# usage:
#########
set help = "\nUsage: hotspot [--help] [--datadir="$datadir"] [--tmpdir="$basetmpdir"/] <input-tags> <outdir>\n"
set help = "$help\t<input-tags> is a per-base bed|starch file with the number of cleavages in the 5th column, excluding bases with no cuts.\n"
set help = "$help\t<outdir> is an output directory where to place badspot calls and minimally (z-score) thresholded hotspot calls.\n"
set help = "$help\tThe [--datadir] should contain the following 3+ column starch files:\n"
set help = "$help\t    blacklist.starch (or .bed)       : Pre-determined troublesome regions to ignore.  A warning is issued if not found.\n"
set help = "$help\t    uniquely-mapping.starch (or .bed): All uniquely-mapping regions, merged.\n"
set help = "$help\tThe [--tmpdir] overrides environmental "\$"TMPDIR which overrides the default "$basetmpdir"/ working temp directory.\n"


#################
# input checking
#################
if ( $#argv == 0 ) then
  printf "$help"
  exit -1
endif

#------------
# want help?
#------------
foreach argc (`seq 1 $#argv`)
  if ( "$argv[$argc]" == "--help" ) then
    printf "$help"
    exit 0
  endif
end

#----------------------
# overriding defaults?
#----------------------
@ nargs = $#argv
@ arglimit = `echo "$#argv-2" | bc -q`
foreach argc (`seq 1 $arglimit`)
  if ( "$argv[$argc]" =~ "--datadir=*" ) then
    set datadir = `echo "$argv[$argc]" | cut -f2 -d'='`
    if ( ! $%datadir ) then
      printf "No value given for --datadir=value\n"
      exit -1
    else if ( ! -s $datadir || ! -d $datadir || ! -r $datadir ) then
      printf "Data directory does not exist or cannot be read: %s\n" $datadir
      exit -1
    endif
    @ nargs--
  else if ( "$argv[$argc]" =~ "--tmpdir=*" ) then
    set tmpdir = `echo "$argv[$argc]" | cut -f2 -d'='`
    if ( ! $%tmpdir ) then
      printf "No value given for --tmpdir=value\n"
      exit -1
    endif
    @ nargs--
    @ defaulttmpdir = 0
  endif
end

#-----------------
# required inputs
#-----------------
@ nreqargs = 2
if ( $nargs != $nreqargs ) then
  printf "$help"
  printf "Wrong number of arguments.  Expect %s mandatory arguments.\n" $nreqargs
  exit -1
endif

@ arglimit++
set tags   = $argv[$arglimit]
@ arglimit++
set outdir = $argv[$arglimit]

if ( ! -s $tags ) then
  printf "Unable to find <input-tags>: %s\n" $tags
  exit -1
else if ( -e $outdir && ! -d $outdir ) then
  printf "<outdir> was set to an existing NON-directory at %s\n" $outdir
  exit -1
endif
mkdir -p $outdir

if ( $defaulttmpdir > 0 ) then
  rm -rf $tmpdir
endif
mkdir -p $tmpdir


####################
# default hardcodes
####################
@ bgnd_half_range      = 25000 # background half-window size
@ fg_half_win          = 125   # foreground half-window size
@ min_hotspot_size     = 10    # smallest hotspot allowed
@ within_range         = 75    # if two hotspots are within 2*$within_range, merge them
@ min_zscore           = 2     # hotspot zscore must be >= $min_zscore
@ n_passes             = 2     # number of hotspot passes, removing all tags from previous pass(es)
@ badspot_mintags      = 5     # consider badspot if local window has at least this many tags
set badspot_threshold  = 0.8   # is a badspot if this proportion of tags occurs in the small window
@ bspot_small_half_win = 25    # badspot foreground half-window size
@ bspot_big_half_win   = 125   # badspot background half-window size
set black_outs         = $datadir/blacklist.starch
set uniqs_maps         = $datadir/uniquely-mapping.starch
set error_log          = $outdir/error.log
set outfile            = $outdir/hotspot.zscore-$min_zscore.starch

if ( $fg_half_win >= $bgnd_half_range ) then
  printf "Your <background-half-range> should be larger than %s\n" $fg_half_win
  exit -1
else if ( ! -s $uniqs_maps ) then
  set uniqs_maps = $datadir/uniquely-mapping.bed
  if ( ! -s $uniqs_maps ) then
    printf "Unable to find <uniq-mapability>: %s in %s\n" $uniqs_maps:t $datadir
    exit -1
  endif
endif

if ( ! -s $black_outs ) then
  set black_outs = $datadir/blacklist.bed
  if ( ! -s $black_outs ) then
    printf "Warning: Expect a file of known suspect regions named %s in %s\n" $black_outs:t $datadir >> $error_log
    printf "\tContinuing anyway...\n" >> $error_log
    set black_outs = ""
  endif
endif


#####################
# determine badspots
#####################
rm -f $tmpdir/.empty
touch $tmpdir/.empty
set badstuff       = ($tmpdir/.empty $black_outs)
set bspt_small_win = $tmpdir/badspot.win-small
set bspt_small_rng = $tmpdir/badspot.win-small-range
set bspt_big_win   = $tmpdir/badspot.win-large
set bspt_big_rng   = $tmpdir/badspot.win-large-range
set badspots       = $outdir/badspots.bed

rm -f $bspt_small_win
rm -f $bspt_big_win
rm -f $bspt_small_rng
rm -f $bspt_big_rng
mkfifo $bspt_small_win
mkfifo $bspt_big_win
mkfifo $bspt_small_rng
mkfifo $bspt_big_rng

# remove known problematic regions from $tags
# use awk '$5 != 0' in case input includes all zeroes in the genome
# use $bspot_small_half_win when offsetting from tags in both small/big windows
#   -> go half a window out from a tag and pad +/- small_win and big_win -> captures the original tag in both cases
(bedops -n 1 $tags $badstuff \
  | awk '$5 != 0' \
  | bedops -u --range $bspot_small_half_win":"$bspot_small_half_win - \
 >! $bspt_small_rng) &

(bedops -n 1 $tags $badstuff \
  | awk '$5 != 0' \
  | bedops -u --range $bspot_small_half_win":"$bspot_small_half_win - \
 >! $bspt_big_rng) &

(bedops -n 1 $tags $badstuff \
  | awk '$5 != 0' \
  | bedmap --faster --prec 0 --delim "\t" --range $bspot_small_half_win --echo-map-range --sum $bspt_small_rng - \
 >! $bspt_small_win) &

(bedops -n 1 $tags $badstuff \
  | awk '$5 != 0' \
  | bedmap --faster --prec 0 --range $bspot_big_half_win --sum $bspt_big_rng - \
 >! $bspt_big_win) &

paste $bspt_small_win $bspt_big_win \
  | awk -v m=$badspot_mintags -v t=$badspot_threshold '($(NF-1) >= m && $(NF-1)/$(NF) >= t)' \
  | sort-bed --max-mem 2G - \
  | bedops -m - \
 >! $badspots

wait

set badstuff = ($badstuff $badspots)

# clean up pipes
rm -f $bspt_small_win
rm -f $bspt_big_win
rm -f $bspt_small_rng
rm -f $bspt_big_rng


#########################################
# run the hotspot algorithm on the input
#########################################
set accumulated = ()
foreach iteration (`seq 1 $n_passes`)
  set bgnd       = $tmpdir/bkgd
  set bgnd_umaps = $tmpdir/umaps

  rm -f $bgnd
  rm -f $bgnd_umaps

  mkfifo $bgnd
  mkfifo $bgnd_umaps

  (bedops -n 1 $tags $badstuff $accumulated \
    | awk '$5 != 0' \
    | bedops -u --range $bgnd_half_range - \
    | bedmap --faster --bases-uniq - $uniqs_maps \
   >! $bgnd_umaps) &

  (bedops -n 1 $tags $badstuff $accumulated \
    | awk '$5 != 0' \
    | bedmap --faster --prec 0 --sum --range $bgnd_half_range - \
   >! $bgnd) &

  bedops -n 1 $tags $badstuff $accumulated \
    | awk '$5 != 0' \
    | cut -f1-5 \
    | bedmap --faster --echo --prec 0 --sum --range $fg_half_win --delim "\t" - \
    | cut -f1-3,6 \
    | paste - $bgnd $bgnd_umaps \
    | awk -v sz=$fg_half_win \
        'BEGIN {OFS="\t"; fsz=2*sz} ; { \
          if ( $6 > fsz ) { \
            p=fsz/$6; \
            expect=$5*p; \
            deviant=sqrt(expect*(1-p)); \
            fgnd_sum=$4; \
            zscore=(fgnd_sum-expect)/deviant; \
            print $1, $2, $3, ".", zscore; \
          } else { \
            print $1, $2, $3, ".", -1; \
          } \
        }' \
    | bedmap --faster --echo --min --range $fg_half_win --delim "\t" - \
    | awk -v minz=$min_zscore \
        'BEGIN {OFS="\t"; lstchrom=""; chrom=""; start=1; end=0; zsc=-1; in_hs=0} ; { \
          if ( $(NF) >= minz ) { \
            if ( in_hs && $1 == lstchrom ) { \
              end = $3; \
              if ( $5 > zsc ) { zsc = $5; } \
            } else { \
              if ( in_hs ) { \
                print chrom, start, end, "i", zsc; \
              } \
              lstchrom = $1; \
              chrom = $1; \
              start = $2; \
              end = $3; \
              zsc = $5; \
            } \
            in_hs = 1; \
          } else { \
            if ( in_hs ) { \
              print chrom, start, end, "i", zsc; \
              start=1; end=0; zsc=-1; chrom=$1; lstchrom=""; \
            } \
            in_hs = 0; \
          } \
        } END { \
          if ( in_hs ) { print chrom, start, end, "i", zsc; } \
        }' \
    | awk -v mn=$min_hotspot_size '$3-$2 >= mn' \
   >! $tmpdir/iteration.$iteration

  wait

  set accumulated = ($accumulated $tmpdir/iteration.$iteration)

  # clean up pipes
  rm -f $bgnd
  rm -f $bgnd_umaps
end

# create file output
set merged_hotspots = $tmpdir/merged.hotspots
rm -f $merged_hotspots
mkfifo $merged_hotspots

(bedops -m --range $within_range $accumulated \
  | bedops -u --range -$within_range - \
 >! $merged_hotspots) &

bedops -u $accumulated \
  | bedmap --echo --max --delim "\t" $merged_hotspots - \
  | awk 'BEGIN {OFS="\t"} ; { print $1, $2, $3, "id-"NR, $4 }' \
  | starch - \
 >! $outfile

wait

if ( $defaulttmpdir > 0 ) then
  rm -rf $tmpdir
endif

exit 0


##############################################################################################################################
# changes relative to hotspot v4 at http://www.uwencode.org/proj/hotspot/
#  - backgrounds are symmetrical about every location with a tag rather than pre-chunked in fixed-size bins across the genome
#       this was the original intention of hotspot
#  - badspots are clipped to size using bedmap --echo-map-range
#  - uses a single foreground window of size $fg_half_win instead of choosing the best result over 3 local window sizes
#       it was a group decision to choose one, and Bob Thurman concurred with this approach
#       the single size is the middle of the 3 used in the original hotspot
##############################################################################################################################
