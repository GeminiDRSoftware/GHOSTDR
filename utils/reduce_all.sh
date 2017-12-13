#!/bin/bash -e
set -o pipefail

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case.
COREDIR=$PWD

CHECK=false

BINNING='2x4'
SEEING=0.5

# Now we define the context. For Quality Assessment, use '--qa', for Quick Look
# use '--ql', for Science Quality, leave blank
QUALITY=''

###### NOTE THAT THE SINGLE SEEING IS BEING USED HERE INSTEAD OF BOTH ########

#Now you have the option of pointing to different directories for each file type.
#All in the same place is the way forward.
BIASDIR=$COREDIR
DARKDIR=$COREDIR
FLATDIR=$COREDIR
ARCDIR=$COREDIR
OBJDIR=$COREDIR

#This line is for cleaning up your directory of all the stuff the reduction creates.
#rm *stack* *forStack* adc* *.log *.list tmp* *_dark.fits *_bias.fits GHOST* *_arc.fits *_flat.fits *slit* *slit* *darkCorrected*

allow_inspection() {
    if $CHECK; then
        echo 'You can now check the reduction at this step.'
        read -p "Press any key to continue... " -n1 -s
    fi
}

add_to_calib_mgr() {
    grep 'Calibration stored as' - | awk '{print $4}' | {
        read calib
        caldb add -v $calib
    }
}

rm -rf calibrations .reducecache  # prepare_data.py outputs not needed because
caldb init -v -w  # we're using the local calibration manager instead

echo 'Splitting MEFs'
#typewalk --tags BUNDLE -n -o bundle
#reduce --drpkg ghostdr @bundle

echo 'Doing slits now'

typewalk --tags GHOST SLITV BIAS UNPREPARED --dir $BIASDIR/ -n -o bias.list
reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)
allow_inspection

typewalk --tags GHOST SLITV DARK UNPREPARED --dir $DARKDIR/ -n -o dark.list
reduce --drpkg ghostdr $QUALITY @dark.list 2>&1 | tee >(add_to_calib_mgr)
allow_inspection

for MODE in HIGH STD; do
    typewalk --tags GHOST SLITV FLAT UNPREPARED $MODE --dir $FLATDIR/ -n -o flat.list
    reduce --drpkg ghostdr $QUALITY @flat.list 2>&1 | tee >(add_to_calib_mgr)
    allow_inspection
    
    while read ARC <&3; do
        echo Reducing slit arc $ARC
        reduce --drpkg ghostdr $QUALITY $ARC 2>&1 | tee >(add_to_calib_mgr)
        allow_inspection
    done 3< <( typewalk --tags GHOST SLITV ARC UNPREPARED $MODE --dir $ARCDIR/ -n -o /dev/fd/2 2>&1 1>/dev/null | grep -v '^#' )

    
    while read STAND <&3; do
        echo Reducing slit standard $STAND
        reduce --drpkg ghostdr $QUALITY $STAND 2>&1 | tee >(add_to_calib_mgr)
        allow_inspection
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE UNPREPARED $MODE --dir $OBJDIR/ --filemask "standard.*\.(fits|FITS)" -n -o /dev/fd/2 2>&1 1>/dev/null \
        | grep -v '^#'
    )

    while read OBJECT <&3; do
        echo Reducing slit object $OBJECT
        reduce --drpkg ghostdr $QUALITY $OBJECT 2>&1 | tee >(add_to_calib_mgr)
        allow_inspection
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE UNPREPARED $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*\.(fits|FITS)" -n -o /dev/fd/2 2>&1 1>/dev/null \
        | grep -v '^#'
    )
done

echo 'Now the spectrograph data'

for CAM in RED BLUE; do
    echo "Doing $CAM images now"

    typewalk --tags GHOST BIAS UNPREPARED 1x1 $CAM --dir $BIASDIR/ -n -o bias.list
    reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)
    allow_inspection
    
    if [ $BINNING != '1x1' ]; then
        typewalk --tags GHOST BIAS UNPREPARED $BINNING $CAM --dir $BIASDIR/ -n -o bias.list
        reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)
        allow_inspection
    fi
    
    typewalk --tags GHOST DARK UNPREPARED $CAM --dir $DARKDIR/ -n -o dark.list
    reduce --drpkg ghostdr $QUALITY @dark.list 2>&1 | tee >(add_to_calib_mgr)
    allow_inspection
    
    for MODE in HIGH STD; do
        typewalk --tags GHOST FLAT UNPREPARED $CAM $MODE --dir $FLATDIR/ -n -o flat.list
        reduce --drpkg ghostdr $QUALITY @flat.list 2>&1 | tee >(add_to_calib_mgr)
        allow_inspection

        while read ARC <&3; do
            echo Reducing arc $ARC
            reduce --drpkg ghostdr $QUALITY $ARC 2>&1 | tee >(add_to_calib_mgr)
            allow_inspection
        done 3< <( typewalk --tags GHOST ARC UNPREPARED $CAM $MODE --dir $ARCDIR/ -n -o /dev/fd/2 2>&1 1>/dev/null | grep -v '^#' )

        while read STAND <&3; do
            echo Reducing standard $STAND
            reduce --drpkg ghostdr $QUALITY $STAND
            allow_inspection
        done 3< <(
            typewalk --tags GHOST UNPREPARED $BINNING $CAM $MODE --dir $OBJDIR/ --filemask "standard.*\.(fits|FITS)" -n -o /dev/fd/2 2>&1 1>/dev/null \
            | grep -v '^#'
        )

        while read OBJECT <&3; do
            echo Reducing object $OBJECT
            reduce --drpkg ghostdr $QUALITY $OBJECT
            allow_inspection
        done 3< <(
            typewalk --tags GHOST UNPREPARED $BINNING $CAM $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*\.(fits|FITS)" -n -o /dev/fd/2 2>&1 1>/dev/null \
            | grep -v '^#'
        )
    done
done
