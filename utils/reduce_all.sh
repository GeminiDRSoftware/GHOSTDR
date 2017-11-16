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

BINNING='1x1'
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

add_to_calib_mgr() {
    grep 'Calibration stored as' - | awk '{print $4}' | {
        read calib
        reduce_db.py add -v $calib
    }
}

rm -rf calibrations .reducecache  # prepare_data.py outputs not needed because
reduce_db.py init -v -w  # we're using the local calibration manager instead

echo 'Doing slits now'

typewalk --tags GHOST SLITV BIAS UNPREPARED --dir $BIASDIR/ -n -o bias.list
reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)

if $CHECK; then
    echo 'You can now check the reduction at this step.'
    read -p "Press any key to continue... " -n1 -s
fi    

typewalk --tags GHOST SLITV DARK UNPREPARED --dir $DARKDIR/ -n -o dark.list
reduce --drpkg ghostdr $QUALITY @dark.list 2>&1 | tee >(add_to_calib_mgr)

if $CHECK; then
    echo 'You can now check the reduction at this step.'
    read -p "Press any key to continue... " -n1 -s
fi

for MODE in HIGH STD; do
    typewalk --tags GHOST SLITV FLAT UNPREPARED $MODE --dir $FLATDIR/ -n -o flat.list
    reduce --drpkg ghostdr $QUALITY @flat.list 2>&1 | tee >(add_to_calib_mgr)

    if $CHECK; then
        echo 'You can now check the reduction at this step.'
        read -p "Press any key to continue... " -n1 -s
    fi
    
    typewalk --tags GHOST SLITV ARC UNPREPARED $MODE --dir $ARCDIR/ -n -o arc.list
    reduce --drpkg ghostdr $QUALITY @arc.list 2>&1 | tee >(add_to_calib_mgr)

    if $CHECK; then
        echo 'You can now check the reduction at this step.'
        read -p "Press any key to continue... " -n1 -s
    fi
    
    while read object <&3; do
        echo Reducing $object
        reduce --drpkg ghostdr $QUALITY $OBJDIR/$object 2>&1 | tee >(add_to_calib_mgr)
        if $CHECK; then
            echo 'You can now check the reduction at this step.'
            read -p "Press any key to continue... " -n1 -s
        fi
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE UNPREPARED $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*\.(fits|FITS)" -n \
        | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" \
        | grep '^\s' \
        | awk '{print $1}'
    )
done

echo 'Now the spectrograph data'

for CAM in RED BLUE; do
    echo "Doing $CAM images now"

    typewalk --tags GHOST BIAS UNPREPARED 1x1 $CAM --dir $BIASDIR/ -n -o bias.list
    reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)
    
    if [ $BINNING != '1x1' ]; then
        typewalk --tags GHOST BIAS UNPREPARED $BINNING $CAM --dir $BIASDIR/ -n -o bias.list
        reduce --drpkg ghostdr $QUALITY @bias.list 2>&1 | tee >(add_to_calib_mgr)
    fi

    if $CHECK; then
        echo 'You can now check the reduction at this step.'
        read -p "Press any key to continue... " -n1 -s
    fi
    
    typewalk --tags GHOST DARK UNPREPARED $CAM --dir $DARKDIR/ -n -o dark.list
    reduce --drpkg ghostdr $QUALITY @dark.list 2>&1 | tee >(add_to_calib_mgr)
    
    if $CHECK; then
        echo 'You can now check the reduction at this step.'
        read -p "Press any key to continue... " -n1 -s
    fi
    
    for MODE in HIGH STD; do
        typewalk --tags GHOST FLAT UNPREPARED $CAM $MODE --dir $FLATDIR/ -n -o flat.list
        reduce --drpkg ghostdr $QUALITY @flat.list 2>&1 | tee >(add_to_calib_mgr)

        if $CHECK; then
            echo 'You can now check the reduction at this step.'
            read -p "Press any key to continue... " -n1 -s
        fi

        typewalk --tags GHOST ARC UNPREPARED $CAM $MODE --dir $ARCDIR/ -n -o arc.list
        reduce --drpkg ghostdr $QUALITY @arc.list  2>&1 | tee >(add_to_calib_mgr)

        if $CHECK; then
            echo 'You can now check the reduction at this step.'
            read -p "Press any key to continue... " -n1 -s
        fi

        while read object <&3; do
            echo Reducing $object
            reduce --drpkg ghostdr $QUALITY $OBJDIR/$object
            if $CHECK; then
                echo 'You can now check the reduction at this step.'
                read -p "Press any key to continue... " -n1 -s
            fi
        done 3< <(
            typewalk --tags GHOST UNPREPARED $BINNING $CAM $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*\.(fits|FITS)" -n \
            | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" \
            | grep '^\s' \
            | awk '{print $1}'
        )
    done
done
