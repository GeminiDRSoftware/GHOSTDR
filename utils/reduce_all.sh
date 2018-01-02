#!/bin/bash -e
# Script to reduce all GHOST data using the recipe system.
# This must be ran from within where the GHOST data are kept

set -o pipefail
exec 6<&0

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'rm -rf /tmp/$$.mark; kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

BINNING=1x1
SEEING=0.5  # Note that the single seeing is being used here instead of both
QUALITY=  # Quality Assessment = --qa, Quick Look = --ql, Science Quality = leave blank
CHECK=true  # pause (true) or not (false) after each 'reduce' call to inspect results
LINGER="${1:-0}"  # how many secs to pause between 'reduce' calls; 1st script arg or 0 default

allow_inspection() {
	builtin echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s -u6
}

# do things after 'reduce': ingest calibs into calmgr, allow inspections, delete debris files
postp() {
	read calib && { [ -f "$calib" ] && caldb add -v $calib; }
	$CHECK && allow_inspection
	# TODO: This is a useful line, but it gets rid of all fits files created after
	# reducing anything. These include the desired pipeline products too, so those should remain or
	# be copied into a 'reduced' directory.
	#find -maxdepth 1 -newer /tmp/$$.mark -type f -name "*.fits" -exec rm -rvf '{}' + 2>/dev/null || true
	$CHECK || sleep $LINGER
}

# do things before 'reduce': save timestamp (used to delete debris files), print a visually obvious banner
prep() {
	[[ "$@" =~ BUNDLE ]] || touch /tmp/$$.mark
	cat <<-HERE  # print the banner/header
		
		
		
		
		
		
		
		    > --------------------------------------------------------------------------
		    >
		    >   $1
		    >
		    > --------------------------------------------------------------------------
		
	HERE
}

# generate the list of files in 'list'; also sets 'msg'
mklist() {
	msg="$1"; shift
	list=`typewalk --tags GHOST UNPREPARED $@ -n -o /dev/stderr 2>&1 1>/dev/null | grep -v '^#'`
	[ -n "$list" ]
}

# perform the reduction; return any calibrator produced
doreduce() {
	reduce --drpkg ghostdr $QUALITY "$@" 2>&1 | tee /dev/tty | { grep 'Calibration stored as' || true; } | awk '{print $4}'
}

# process all matching files together (as a whole / in a batch)
reduce_list() {
	mklist "$@" || return 0
	prep "$@"
	doreduce @/dev/stdin <<<"$list" | postp
}

# process each matching file individually (by itself)
reduce_each() {
	mklist "$@" || return 0
	while read THING; do
		[[ "${THING}" =~ .*/(.*) ]] && prep "$msg: ${BASH_REMATCH[1]}"
		doreduce $THING | postp
	done <<<"$list"
}

rm -rf calibrations .reducecache reduce.log  # start with a fresh local cache and logfile
caldb init -v -w  # start with a fresh local calibration manager

reduce_list "Splitting MEFs" BUNDLE  # no need to comment out: noop's on -split simulator outputs
for CAM in SLITV RED BLUE; do
	# process biases (populate a hash with keys for each necessary binning mode, and run 'reduce' for each)
	unset bins; declare -A bins  # 'bins' is the hash
	if [ $CAM = SLITV ]; then bins[2x2]=; else bins[1x1]=; bins[$BINNING]=; fi  # populate
	for bin in "${!bins[@]}"; do reduce_list "Reducing $CAM biases" $CAM BIAS $bin; done  # iterate

	# process everything else
	bin=$BINNING; [ $CAM = SLITV ] && bin=  # binning modes for objects and standards
	reduce_list "Reducing $CAM darks" $CAM DARK
	for MODE in HIGH STD; do
		reduce_list "Reducing $CAM $MODE flats" $CAM $MODE FLAT
		reduce_each "Reducing $CAM $MODE arc" $CAM $MODE ARC
		reduce_each "Reducing $CAM $MODE standard" $CAM $MODE $bin --filemask "standard.*\.(fits|FITS)" 
		reduce_each "Reducing $CAM $MODE object" $CAM $MODE $bin --filemask "obj.*$SEEING.*\.(fits|FITS)" 
	done
done
