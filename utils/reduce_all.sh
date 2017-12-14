#!/bin/bash -e
set -o pipefail
linger="$1"

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'rm -rf /tmp/$$.fifo /tmp/$$.mark; kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept

BINNING='2x4'

###### NOTE THAT THE SINGLE SEEING IS BEING USED HERE INSTEAD OF BOTH ########
SEEING=0.5

# Now we define the context. For Quality Assessment, use '--qa', for Quick Look
# use '--ql', for Science Quality, leave blank
QUALITY=''

mkfifo /tmp/$$.fifo  # IPC between main script and tee'd subprocesses

# change 'false' to 'true' below if you want to inspect outputs between 'reduce' invocations
allow_inspection() {
	if false; then
		builtin echo 'You can now check the reduction at this step.'
		read -p "Press any key to continue... " -n1 -s
	fi
}

# "post-process": call immediately following 'reduce' (must also tee 'reduce' output to
# 'send_cal_to_postp'); this waits on tee'd processes, ingests calibs into calmgr, allows
# inspections, and deletes debris files created by 'reduce'
postp() {
	read calib </tmp/$$.fifo || true
	[ -n "$calib" ] && caldb add -v $calib
	allow_inspection
	find -maxdepth 1 -newer /tmp/$$.mark -type f -name "*.fits" -exec rm -rvf '{}' + 2>/dev/null || true
}

# strip calibrator path (if any) from 'reduce' output and send to 'postp' via fifo
send_cal_to_postp() {
	grep 'Calibration stored as' - | awk '{print $4}' >/tmp/$$.fifo
}

# "pre-process": call before 'reduce'; records timestamp (used later to delete debris
# files), prints large, more obvious header to separate 'reduce' invocations from one another
# (relies on tab indentation so entire script has been converted)
prep() {
	sleep ${linger:=0}
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
	list=`typewalk --tags GHOST UNPREPARED $@ -n -o /dev/fd/2 2>&1 1>/dev/null | grep -v '^#' || true`
	[ -n "$list" ]
}

# process all matching files together (as a whole / in a batch)
reduce_list() {
	mklist "$@" || return 0
	prep "$@"
	reduce --drpkg ghostdr $QUALITY @/dev/fd/0 2>&1 <<<"$list" | tee >(send_cal_to_postp)
	postp
}

# process each matching file individually (by itself)
reduce_each() {
	mklist "$@" || return 0
	while read THING <&3; do
		[[ "${THING}" =~ .*/(.*) ]] && prep "$msg: ${BASH_REMATCH[1]}"
		reduce --drpkg ghostdr $QUALITY $THING 2>&1 | tee >(send_cal_to_postp)
		postp
	done 3<<<"$list"
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
