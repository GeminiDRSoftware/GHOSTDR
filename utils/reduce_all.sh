#!/usr/bin/bash

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

#$DIR = `pwd`

#rm -Rf calibrations/
#
# echo 'Doing slits now'
# typewalk --types GHOST_SLITV_BIAS --dir biases/ -o bias.list
# reduce @bias.list
# typewalk --types GHOST_SLITV_DARK --dir darks/ -o dark.list
# reduce @dark.list --override_cal processed_bias:`ls calibrations/storedcals/bias*SLIT*.fits`

# #FLATS
# typewalk --types GHOST_SLITV_FLAT GHOST_SLITV_HIGH --dir flats/high/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*SLIT*.fits` processed_dark:`ls calibrations/storedcals/dark*SLIT*.fits`
# typewalk --types GHOST_SLITV_FLAT --dir flats/std/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*SLIT*.fits` processed_dark:`ls calibrations/storedcals/dark*SLIT*.fits`

# #ARCS
# typewalk --types GHOST_SLITV_ARC GHOST_SLITV_HIGH --dir arcs/ -o arc.list
# reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*SLIT*.fits` processed_dark:`ls calibrations/storedcals/dark*SLIT*.fits` processed_slitflat:`ls calibrations/storedcals/flat*high*SLIT*.fits` 
# typewalk --types GHOST_SLITV_ARC GHOST_SLITV_STD --dir arcs/ -o arc.list
# reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*SLIT*.fits` processed_dark:`ls calibrations/storedcals/dark*SLIT*.fits` processed_slitflat:`ls calibrations/storedcals/flat*std*SLIT*.fits`



# echo 'Now the spectrograph data'
# echo 'Starting with the blue arm'
# typewalk --types GHOST_BIAS GHOST_BLUE --dir biases/ -o bias.list
# reduce @bias.list
# typewalk --types GHOST_DARK GHOST_BLUE --dir darks/ -o dark.list
# reduce @dark.list --override_cal processed_bias:`ls calibrations/storedcals/bias_*_blue_*.fits`

# #FLATS
# typewalk --types GHOST_FLAT GHOST_BLUE GHOST_HIGH --dir flats/high/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*blue_*.fits` processed_dark:`ls calibrations/storedcals/dark*blue*.fits`
# typewalk --types GHOST_FLAT GHOST_BLUE GHOST_STD --dir flats/std/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*blue_*.fits` processed_dark:`ls calibrations/storedcals/dark*blue*.fits`

#ARCS
typewalk --types GHOST_ARC GHOST_BLUE GHOST_HIGH --dir arcs/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*blue_*.fits` processed_dark:`ls calibrations/storedcals/dark*blue*.fits` processed_slitflat:`ls calibrations/storedcals/flat*high*SLIT*.fits` processed_slit:`ls calibrations/storedcals/arc*high*SLIT*.fits` processed_polyfit:`ls calibrations/storedcals/*blue*high*xmod*.fits`
typewalk --types GHOST_ARC GHOST_BLUE GHOST_STD --dir arcs/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*blue_*.fits` processed_dark:`ls calibrations/storedcals/dark*blue*.fits` processed_slitflat:`ls calibrations/storedcals/flat*high*SLIT*.fits` processed_slit:`ls calibrations/storedcals/arc*std*SLIT*.fits`  processed_polyfit:`ls calibrations/storedcals/*blue*std*xmod*.fits`




# echo 'Doing red now'
# typewalk --types GHOST_BIAS GHOST_RED --dir biases/ -o bias.list
# reduce @bias.list
# typewalk --types GHOST_DARK GHOST_RED --dir darks/ -o dark.list
# reduce @dark.list --override_cal processed_bias:`ls calibrations/storedcals/bias_*_red_*.fits`

# #FLATS
# typewalk --types GHOST_FLAT GHOST_RED GHOST_HIGH --dir flats/high/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*red_*.fits` processed_dark:`ls calibrations/storedcals/dark*red*.fits`
# typewalk --types GHOST_FLAT GHOST_RED GHOST_STD --dir flats/std/ -o flat.list
# reduce @flat.list --override_cal processed_bias:`ls calibrations/storedcals/bias*red_*.fits` processed_dark:`ls calibrations/storedcals/dark*red*.fits`

#ARCS
typewalk --types GHOST_ARC GHOST_RED GHOST_HIGH --dir arcs/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*red_*.fits` processed_dark:`ls calibrations/storedcals/dark*red*.fits` processed_slitflat:`ls calibrations/storedcals/flat*high*SLIT*.fits` processed_slit:`ls calibrations/storedcals/arc*high*SLIT*.fits` processed_polyfit:`ls calibrations/storedcals/*red*high*xmod*.fits`
typewalk --types GHOST_ARC GHOST_RED GHOST_STD --dir arcs/ -o arc.list
reduce @arc.list --override_cal processed_bias:`ls calibrations/storedcals/bias*red_*.fits` processed_dark:`ls calibrations/storedcals/dark*red*.fits` processed_slitflat:`ls calibrations/storedcals/flat*std*SLIT*.fits` processed_slit:`ls calibrations/storedcals/arc*std*SLIT*.fits`  processed_polyfit:`ls calibrations/storedcals/*red*std*xmod*.fits`



