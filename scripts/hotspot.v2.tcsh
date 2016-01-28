#!/bin/tcsh -ef
# author : Shane Neph
# date   : 2016
# proj   : Reimplement hotspot using BEDOPS; see bottom of file for improvements over original hotspot
#            - includes badspot calculations, and follow-on thresholded z-score hotspot calculations
#            - may be called on real or simulated data.


#########
# usage:
#########
set help = "\nUsage: hotspot [--help] [--datadir=../data/hg19] [--tmpdir=/tmp/] <input-tags> <outdir>\n"
set help = "$help\t<input-tags> is a per-base bed|starch file with the number of 5' cuts in the 5th column, excluding bases with 0 cuts.\n"
set help = "$help\t<outdir> is an output directory where to place minimally (z-score) thresholded hotspot calls and badspot calls.\n"
set help = "$help\t[--datadir] should contain the following 3+ column starch files:\n"
set help = "$help\t    blacklist.starch        : Pre-determined troublesome regions to ignore.  A warning is issued if not found.\n"
set help = "$help\t    uniquely-mapping.starch : All uniquely-mapping regions, merged.\n"
set help = "$help\t[--tmpdir] overrides environmental "\$"TMPDIR which overrides the default /tmp/ working temp directory.\n"


###########
# defaults
###########
set db      = hg19
set datadir = ../data/$db
set tmpdir  = /tmp/`whoami`/hotspot/$$
if ( $?TMPDIR ) then
  # use the environmental variable
  set tmpdir = $TMPDIR"/$$"
endif


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
foreach argc (`seq 1 $#argv`)
  if ( "$argv[$argc]" =~ "--datadir=*" ) then
    set datadir = `echo "$argv[$argc]" | cut -f2 -d'='`
    if ( ! $%datadir ) then
      printf "No value given for --datadir=value\n"
      exit -1
    else if ( ! -s $datadir ) then
      printf "Data diretory does not exist: %s\n" $datadir
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

set tags   = $1
set outdir = $2

if ( ! -s $tags ) then
  printf "Unable to find <input-tags>: %s\n" $tags
  exit -1
else if ( -e $outdir && ! -d $outdir ) then
  printf "<outdir> was set to an existing NON-directory at %s\n" $outdir
  exit -1
endif
mkdir -p $outdir


############
# hardcodes
############
@ bgnd_half_range      = 25000 # background half-window size
@ fg_half_win          = 125   # foreground half-window size
@ min_hotspot_size     = 10    # smallest hotspot allowed
@ within_range         = 75    # if two hotspots are within 2*$within_range, merge them
@ min_zscore           = 2     # hotspot zscore must be >= $min_zscore
@ n_passes             = 2     # number of hotspot passes, removing all tags from previous pass(es)
@ badspot_mintags      = 5     # consider badspot if local window has at least this many tags
@ badspot_threshold    = 0.8   # is a badspot if this proportion of tags occurs in the small window
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
  printf "Unable to find <uniq-mapability>: %s in %s\n" $uniqs_maps:t $datadir
  exit -1
endif

if ( ! -s $black_outs ) then
  printf "Warning: Expect a file of known suspect regions named %s in %s\n" $black_outs:t $datadir > $errorlog
  printf "\tContinuing anyway...\n" > $errorlog
  set black_outs = ""
endif


############################
# temp directory to do work
############################
rm -rf $tmpdir
mkdir -p $tmpdir


#####################
# determine badspots
#####################
set bspt_small_win = $tmpdir/badspot.win-small
set bspt_big_win   = $tmpdir/badspot.win-large
set badspots       = $tmpdir/badspots.bed

rm -f $bspt_small_win
rm -f $bspt_big_win
mkfifo $bspt_small_win
mkfifo $bspt_big_win

(bedops -u --range $bspot_small_half_win $tags \
  | bedmap --faster --prec 0 --delim "\t" --echo --sum - \
 >! $bspt_small_win) &

(bedops -u --range $bspot_big_half_win $tags \
  | bedmap --faster --prec 0 --sum - \
 >! $bspt_big_win) &

paste $bspt_small_win $bspt_big_win \
  | awk -v m=$badspot_mintags -v t=$badspot_threshold '($(NF-1) >= m && $(NF-1)/$(NF) >= t)' \
  | bedops -m - \
 >! $badspots

wait

# clean up pipes
rm -f $bspt_small_win
rm -f $bspt_big_win


#########################################
# run the hotspot algorithm on the input
#########################################
rm -f $tmpdir/.empty
touch $tmpdir/.empty # for $accumulated in the event that both $black_outs & $badspots dne

set accumulated = ($tmpdir/.empty $black_outs $badspots)
foreach iteration (`seq 1 $n_passes`)
  set bgnd       = $tmpdir/bkgd
  set bgnd_umaps = $tmpdir/umaps

  rm -f $bgnd
  rm -f $bgnd_umaps

  mkfifo $bgnd
  mkfifo $bgnd_umaps

  (bedops -n 1 $tags $accumulated \
    | bedops -u --range $bgnd_half_range - \
    | bedmap --faster --bases-uniq - $uniqs_maps \
   >! $bgnd_umaps) &

  (bedops -n 1 $tags $accumulated \
    | bedmap --faster --prec 0 --sum --range $bgnd_half_range - \
   >! $bgnd) &

  bedops -n 1 $tags $accumulated \
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

rm -rf $tmpdir

exit 0


##############################################################################################################################
# changes relative to hotspot v4 at http://www.uwencode.org/proj/hotspot/
#  - backgrounds are symmetrical about every location with a tag rather than pre-chunked in fixed size bins across the genome
#       this was the original intention of hotspot
#  - uses a single foreground window of size $fg_half_win instead of choosing the best result over 3 local window sizes
#       it was a group decision to choose one, and Bob Thurman concurred in private communication
#       the single size is the middle of the 3 used in the original hotspot
##############################################################################################################################
