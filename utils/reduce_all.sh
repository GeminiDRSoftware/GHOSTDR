#!/bin/sh

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case. 
COREDIR=`pwd`

CALDIR=$COREDIR'/calibrations/storedcals'

BIASDIR=$COREDIR'/biases'
DARKDIR=$COREDIR'/darks'
#Assumes the high and std flats are inside high/ and std/ dirs in FLATDIR
FLATDIR=$COREDIR'/flats'
ARCDIR=$COREDIR'/arcs'


echo 'Doing slits now'
typewalk --types GHOST_SLITV_BIAS --dir $BIASDIR/ -o bias.list
reduce @bias.list
typewalk --types GHOST_SLITV_DARK --dir $DARKDIR/ -o dark.list
reduce @dark.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits`

#FLATS
typewalk --types GHOST_SLITV_FLAT GHOST_SLITV_HIGH --dir $FLATDIR/high/ -o flat.list
reduce @flat.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits`
# TODO: Note that GHOST_SLITV_STD is not a thing, so this needs fixing.
typewalk --types GHOST_SLITV_FLAT --dir $FLATDIR/std/ -o flat.list
reduce @flat.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits`

#ARCS
typewalk --types GHOST_SLITV_ARC GHOST_SLITV_HIGH --dir $ARCDIR/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits` processed_slitflat:`ls $CALDIR/flat*high*SLIT*.fits` 
typewalk --types GHOST_SLITV_ARC GHOST_SLITV_STD --dir $ARCDIR/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls $CALDIR/bias*SLIT*.fits` processed_dark:`ls $CALDIR/dark*SLIT*.fits` processed_slitflat:`ls $CALDIR/flat*std*SLIT*.fits`


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
	typewalk --types GHOST_FLAT GHOST_$CAPCAM GHOST_$CAPMODE --dir $FLATDIR/$j/ -o flat.list
	reduce @flat.list --override_cal processed_bias:`ls $CALDIR/bias*$i*.fits` processed_dark:`ls $CALDIR/dark*$i*.fits`
	typewalk --types GHOST_ARC GHOST_$CAPCAM GHOST_$CAPMODE --dir $ARCDIR/ -o arc.list
	reduce @arc.list --override_cal processed_bias:`ls $CALDIR/bias*$i*.fits` processed_dark:`ls $CALDIR/dark*$i*.fits` processed_slitflat:`ls $CALDIR/flat*$j*SLIT*.fits` processed_slit:`ls $CALDIR/arc*$j*SLIT*.fits` processed_polyfit:`ls $CALDIR/*$i*$j*xmod*.fits`
    done
done





