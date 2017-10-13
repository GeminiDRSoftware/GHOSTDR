#!/bin/bash

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case.
COREDIR=$PWD'/original'

CALDIR=$PWD'/calibrations/storedcals'

BINNING='1x1'
#Now you have the option of pointing to different directories for each file type.
#All in the same place is the way forward.
BIASDIR=$COREDIR
DARKDIR=$COREDIR
FLATDIR=$COREDIR
ARCDIR=$COREDIR
OBJDIR=$COREDIR

#This line is for cleaning up your directory of all the stuff the reduction creates.
#rm *stack* *forStack* adc* *.log *.list tmp* *_dark.fits *_bias.fits GHOST* *_arc.fits *_flat.fits


 echo 'Doing slits now'

typewalk --tags GHOST SLITV BIAS --dir $BIASDIR/ -o bias.list
reduce @bias.list

typewalk --tags GHOST SLITV DARK --dir $DARKDIR/ -o dark.list
reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits`

for mode in high std; do
    CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

    typewalk --tags GHOST SLITV FLAT $CAPMODE --dir $FLATDIR/ -o flat.list
    reduce @flat.list
    
    typewalk --tags GHOST SLITV ARC $CAPMODE --dir $ARCDIR/ -o arc.list
    reduce @arc.list
    
    while read object <&3; do
        echo Reducing $object
        reduce $object
	
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE $CAPMODE --dir $OBJDIR/ --filemask 'obj.*\.(fits|FITS)' \
            -o tmp$$.list >& /dev/null && cat tmp$$.list | grep -v '^#'; rm tmp$$.list
    )
done

echo 'Now the spectrograph data'

for cam in red blue; do
    echo "Doing $cam images now"
    CAPCAM=`echo $cam | tr '[:lower:]' '[:upper:]'`

    typewalk --tags GHOST BIAS $CAPCAM --dir $BIASDIR/ --filemask '.*'$BINNING'.*\.(fits|FITS)' -o bias.list
    reduce @bias.list

    typewalk --tags GHOST DARK $CAPCAM --dir $DARKDIR/ -o dark.list
    reduce @dark.list 

    for mode in high std; do
        CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

        typewalk --tags GHOST FLAT $CAPCAM $CAPMODE --dir $FLATDIR/ -o flat.list
        reduce @flat.list

        typewalk --tags GHOST ARC $CAPCAM $CAPMODE --dir $ARCDIR/ -o arc.list
        reduce @arc.list 

        for seeing in 0.5 1.0; do
            while read object <&3; do
                echo Reducing $object
                reduce $object
		
            done 3< <(
                typewalk --tags GHOST $CAPCAM GHOST_$CAPMODE --dir $OBJDIR/ --filemask ".*$seeing.*'$BINNING'.*\.(fits|FITS)" \
                    -o tmp$$.list >& /dev/null && cat tmp$$.list | grep -v '^#'; rm tmp$$.list
            )
        done
    done
done
