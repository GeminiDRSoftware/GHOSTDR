#!/bin/bash -e
# Script to reduce all GHOST data using the recipe system.
# This must be ran from within where the GHOST data are kept

LC_ALL=C
set -o pipefail
exec 6<&0

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'rm -rf /tmp/$$.mark; kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

BINNING=1x1
SEEING=0.5  # Note that the single seeing is being used here instead of both
QUALITY=  # Quality Assessment = --qa, Quick Look = --ql, Science Quality = leave blank
CHECK=false  # pause (true) or not (false) after each 'reduce' call to inspect results
LINGER="${1:-0}"  # how many secs to pause between 'reduce' calls; 1st script arg or 0 default
DELINT=true  # delete intermediate files (true) or move them into parallel folder (false)
$DELINT || INTERMED=`mktemp -d intermediates.XXXXXXXXXX`

SHORTSPEC=hd200654
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LONGSPEC=$DIR/../simulator/pyghost/data/standards/${SHORTSPEC}.fits

allow_inspection() {
	builtin echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s -u6
}

# do things after 'reduce': ingest calibs into calmgr, allow inspections, delete debris files
postp() {
	read calib && { [ -f "$calib" ] && { echo caldb add -v $calib >>commands && caldb add -v $calib; }; }
	$CHECK && allow_inspection
	[[ "$@" =~ BUNDLE || ( ( "$@" =~ object || "$@" =~ standard ) && ! "$@" =~ SLIT ) ]] || {
		{ find . -maxdepth 1 -newer /tmp/$$.mark -type f -name "*.fits" 2>/dev/null || true; } | {
			if $DELINT; then xargs rm -vf; else xargs -I {} mv -v {} $INTERMED; fi
		}
	}
	$CHECK || sleep $LINGER
}

# do things before 'reduce': save timestamp (used to delete debris files), print a visually obvious banner
prep() {
	touch /tmp/$$.mark
	cat <<-HERE







		    > --------------------------------------------------------------------------
		    >
		    >   $1
		    >
		    > --------------------------------------------------------------------------

	HERE
}

# generate the list of files in 'list'; also sets 'msg'
mklist() {
	flag="$1"; shift
	if [[ $flag == silent ]]; then msg="$1"; shift; else msg=$flag; fi
    if [[ -x /usr/bin/uuidgen ]]; then
        UUID=`/usr/bin/uuidgen`
    else
        UUID=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1 || true`
    fi
	[[ $flag == silent ]] || echo typewalk --adpkg ghost_instruments --tags GHOST UNPREPARED $@ -n -o $UUID >>commands
	list=`typewalk --adpkg ghost_instruments --tags GHOST UNPREPARED $@ -n -o /dev/stderr 2>&1 1>/dev/null | grep -v '^#' || true`
	[ -n "$list" ]
}

# perform the reduction; return any calibrator produced
doreduce() {
	if [[ "$1" == @* ]]; then TARGET=@$UUID; else TARGET=$1; fi
	echo PYSYN_CDBS=. reduce --drpkg ghostdr --adpkg ghost_instruments $STANDARD $QUALITY $TARGET >>commands
	stdbuf -o 0 reduce --drpkg ghostdr $STANDARD $QUALITY "$@" 2>&1 | stdbuf -o 0 grep -v stdbuf | tee /dev/tty \
		| { grep 'Calibration stored as' || true; } | awk '{print $4}'
}

# process all matching files together (as a whole / in a batch)
reduce_list() {
	mklist "$@" || return 0
	prep "$@"
	doreduce @/dev/stdin <<<"$list" | postp "$@"
}

# process each matching file individually (by itself)
reduce_each() {
	mklist silent "$@" || return 0
	while read THING; do
		[[ "${THING}" =~ .*/(.*) ]] && prep "$msg: ${BASH_REMATCH[1]}"
		doreduce $THING | postp "$@"
	done <<<"$list"
	STANDARD=
}

rm -rf calibrations .reducecache reduce.log  # start with a fresh local cache and logfile
caldb init -v -w  # start with a fresh local calibration manager
echo caldb init -v -w >>commands

reduce_list "Splitting MEFs" BUNDLE  # no need to comment out: noop's on -split simulator outputs
for CAM in SLITV BLUE RED; do
	# process biases (populate an array with each necessary binning mode, and run 'reduce' for each)
	
	bins=()  # 'bins' is the array
	if [ $CAM = SLITV ]; then bins+=(2x2); else bins+=($BINNING); [[ "${bins[@]}" =~ 1x1 ]] || bins+=(1x1); fi  # populate
	for BIN in "${bins[@]}"; do reduce_list "Reducing $CAM biases" $CAM BIAS $BIN; done  # iterate
    
    # process everything else
	BIN=$BINNING; [ $CAM = SLITV ] && BIN=  # binning modes for objects and standards
	reduce_list "Reducing $CAM darks" $CAM DARK
	
	echo "Darks Reduced for Cam " + $CAM
	for MODE in HIGH STD; do
		reduce_list "Reducing $CAM $MODE flats" $CAM $MODE FLAT
		reduce_each "Reducing $CAM $MODE arc" $CAM $MODE ARC
		reduce_each "Reducing $CAM $MODE standard" $CAM $MODE $BIN --filemask "standard.*\.(fits|FITS)" 
		STANDARD=`typewalk --adpkg ghost_instruments --tags GHOST $CAM $MODE $BIN --filemask "standard.*${SHORTSPEC}.*wavelengthAdded\.(fits|FITS)" -n -o /dev/stderr 2>&1 1>/dev/null | grep -v '^#' || true`
		STANDARD=${STANDARD:+-p std=$STANDARD std_spec=$LONGSPEC}
		reduce_each "Reducing $CAM $MODE object" $CAM $MODE $BIN --filemask "obj.*$SEEING.*\.(fits|FITS)" 
	done
	echo "Reduction Complete for Cam " + $CAM
done
