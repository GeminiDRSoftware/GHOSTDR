#!/bin/sh

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case. 
COREDIR=$PWD

CALDIR=$COREDIR'/calibrations/storedcals'

#Now you have the option of pointing to different directories for each file type.
#All in the same place is the way forward.
BIASDIR=$COREDIR
DARKDIR=$COREDIR
FLATDIR=$COREDIR
ARCDIR=$COREDIR

#This line is for cleaning up your directory of all the stuff the reduction creates. 
#rm *stack* *forStack* adc* *.log *.list tmp* *_dark.fits *_bias.fits GHOST* *_arc.fits *_flat.fits

   
echo 'Doing slits now'
typewalk --types GHOST_SLITV_BIAS --dir $BIASDIR/ -o bias.list
reduce @bias.list
typewalk --types GHOST_SLITV_DARK --dir $DARKDIR/ -o dark.list
reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits`

for j in high std
do
    CAPMODE=`echo $j | tr '[:lower:]' '[:upper:]'`
    typewalk --types GHOST_SLITV_FLAT GHOST_$CAPMODE --dir $FLATDIR/ -o flat.list
    reduce @flat.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits`
    typewalk --types GHOST_SLITV_ARC GHOST_$CAPMODE --dir $ARCDIR/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits` processed_slitflat:`ls $CALDIR/flat*$j*SLIT*.fits` 
done

echo 'Now the spectrograph data'
for i in blue red 
do
    echo "Doing $i images now"
    CAPCAM=`echo $i | tr '[:lower:]' '[:upper:]'`
    typewalk --types GHOST_BIAS GHOST_$CAPCAM --dir "$BIASDIR" -o bias.list
    reduce @bias.list
    typewalk --types GHOST_DARK GHOST_$CAPCAM --dir $DARKDIR -o dark.list
    reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*$i*.fits`
    for j in high std
    do
	CAPMODE=`echo $j | tr '[:lower:]' '[:upper:]'`
	typewalk --types GHOST_FLAT GHOST_$CAPCAM GHOST_$CAPMODE --dir $FLATDIR/ -o flat.list
	reduce @flat.list --override_cal processed_bias:`ls $CALDIR/bias*$i*.fits` processed_dark:`ls $CALDIR/dark*$i*.fits`
	typewalk --types GHOST_ARC GHOST_$CAPCAM GHOST_$CAPMODE --dir $ARCDIR/ -o arc.list
	reduce @arc.list --override_cal processed_bias:`ls $CALDIR/bias*$i*.fits` processed_dark:`ls $CALDIR/dark*$i*.fits` processed_slitflat:`ls $CALDIR/flat*$j*SLIT*.fits` processed_slit:`ls $CALDIR/arc*$j*SLIT*.fits` processed_xmod:`ls $CALDIR/*$i*$j*xmod*.fits`
    done
done





