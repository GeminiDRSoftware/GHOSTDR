#!/bin/sh

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case.
COREDIR=$PWD'/original'

CALDIR=$PWD'/calibrations/storedcals'

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

typewalk --types GHOST_SLITV_BIAS --dir $BIASDIR/ -o bias.list
reduce @bias.list

typewalk --types GHOST_SLITV_DARK --dir $DARKDIR/ -o dark.list
reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits`

for mode in high std; do
    CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

    typewalk --types GHOST_SLITV_FLAT GHOST_$CAPMODE --dir $FLATDIR/ -o flat.list
    reduce @flat.list --override_cal \
        processed_bias:`ls $CALDIR/bias*SLIT*.fits` \
        processed_dark:`ls $CALDIR/dark*SLIT*.fits`

    typewalk --types GHOST_SLITV_ARC GHOST_$CAPMODE --dir $ARCDIR/ -o arc.list
    reduce @arc.list --override_cal \
        processed_bias:`ls $CALDIR/bias*SLIT*.fits` \
        processed_dark:`ls $CALDIR/dark*SLIT*.fits` \
        processed_slitflat:`ls $CALDIR/flat*$mode*SLIT*.fits`

    while read object <&3; do
        reduce $object --override_cal \
            processed_bias:`ls $CALDIR/bias*SLIT*.fits` \
            processed_dark:`ls $CALDIR/dark*SLIT*.fits` \
            processed_slitflat:`ls $CALDIR/flat*$mode*SLIT*.fits`
    done 3< <(
        typewalk --types GHOST_SLITV_IMAGE GHOST_$CAPMODE --dir $OBJDIR/ --filemask 'obj.*\.(fits|FITS)' \
        | grep '^ ' \
        | awk '{print $1}'
    )
done

echo 'Now the spectrograph data'

for cam in red blue; do
    echo "Doing $cam images now"
    CAPCAM=`echo $cam | tr '[:lower:]' '[:upper:]'`

    typewalk --types GHOST_BIAS GHOST_$CAPCAM --dir $BIASDIR/ --filemask '.*1x1.*\.(fits|FITS)' -o bias.list
    reduce @bias.list

    typewalk --types GHOST_DARK GHOST_$CAPCAM --dir $DARKDIR/ -o dark.list
    reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*$cam*.fits`

    for mode in high std; do
        CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

        typewalk --types GHOST_FLAT GHOST_$CAPCAM GHOST_$CAPMODE --dir $FLATDIR/ -o flat.list
        reduce @flat.list --override_cal \
            processed_bias:`ls $CALDIR/bias*$cam*.fits` \
            processed_dark:`ls $CALDIR/dark*$cam*.fits` \
			processed_slitflat:`ls $CALDIR/flat*$mode*SLIT*.fits` \
            processed_xmod:`ls $CALDIR/*$cam*$mode*xmod*.fits`

        typewalk --types GHOST_ARC GHOST_$CAPCAM GHOST_$CAPMODE --dir $ARCDIR/ -o arc.list
        reduce @arc.list --override_cal \
            processed_bias:`ls $CALDIR/bias*$cam*.fits` \
            processed_dark:`ls $CALDIR/dark*$cam*.fits` \
            processed_slitflat:`ls $CALDIR/flat*$mode*SLIT*.fits` \
            processed_slit:`ls $CALDIR/arc*$mode*SLIT*.fits` \
            processed_xmod:`ls $CALDIR/*$cam*$mode*xmod*.fits`

        while read object <&3; do
            reduce $object --override_cal \
                processed_bias:`ls $CALDIR/bias*$cam*.fits` \
                processed_dark:`ls $CALDIR/dark*$cam*.fits` \
                processed_slitflat:`ls $CALDIR/flat*$mode*SLIT*.fits` \
                processed_slit:`ls $CALDIR/obj*$mode*SLIT*.fits` \
                processed_xmod:`ls $CALDIR/*$cam*$mode*xmod*.fits` \
                processed_flat:`ls $CALDIR/flat*$mode*$cam*.fits`
        done 3< <(
            typewalk --types GHOST_OBJECT GHOST_$CAPCAM GHOST_$CAPMODE --dir $OBJDIR/ --filemask '.*1x1.*\.(fits|FITS)' \
            | grep '^ ' \
            | awk '{print $1}'
        )
    done
done
