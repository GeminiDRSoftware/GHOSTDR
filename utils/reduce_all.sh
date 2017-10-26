#!/bin/bash

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case.
COREDIR=$PWD

CHECK=false

BINNING='1x2'
SEEING=0.5

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


echo 'Doing slits now'

typewalk --tags GHOST SLITV BIAS --dir $BIASDIR/ -n -o bias.list
reduce --drpkg ghostdr @bias.list

if $CHECK
then
    echo 'You can now check the reduction at this step.'
    read -p "Press any key to continue... " -n1 -s
fi    

typewalk --tags GHOST SLITV DARK --dir $DARKDIR/ -n -o dark.list
reduce --drpkg ghostdr @dark.list

if $CHECK
then
    echo 'You can now check the reduction at this step.'
    read -p "Press any key to continue... " -n1 -s
fi

for MODE in HIGH STD; do
    typewalk --tags GHOST SLITV FLAT $MODE --dir $FLATDIR/ -n -o flat.list
    reduce --drpkg ghostdr @flat.list

    if $CHECK
    then
    	echo 'You can now check the reduction at this step.'
    	read -p "Press any key to continue... " -n1 -s
    fi
    
    typewalk --tags GHOST SLITV ARC $MODE --dir $ARCDIR/ -n -o arc.list
    reduce --drpkg ghostdr @arc.list

    if $CHECK
    then
    	echo 'You can now check the reduction at this step.'
    	read -p "Press any key to continue... " -n1 -s
    fi
    
    while read object <&3; do
        echo Reducing $object
        reduce --drpkg ghostdr $OBJDIR/$object
	if $CHECK
	then
	    echo 'You can now check the reduction at this step.'
	    read -p "Press any key to continue... " -n1 -s
	fi
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*\.(fits|FITS)" -n \
            | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" \
            | grep '^\s' \
            | awk '{print $1}'
    )
done

echo 'Now the spectrograph data'

for CAM in RED BLUE; do
    echo "Doing $CAM images now"

    typewalk --tags GHOST BIAS $CAM --dir $BIASDIR/ --filemask '.*1x1.*\.(fits|FITS)' -n -o bias.list
    reduce --drpkg ghostdr @bias.list
    
    if [ $BINNING != '1x1' ]
    then
	typewalk --tags GHOST BIAS $CAM --dir $BIASDIR/ --filemask '.*'$BINNING'.*\.(fits|FITS)' -n -o bias.list
	reduce --drpkg ghostdr @bias.list
    fi

    if $CHECK
    then
	echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s
    fi
    
    typewalk --tags GHOST DARK $CAM --dir $DARKDIR/ -n -o dark.list
    reduce --drpkg ghostdr @dark.list
    
    if $CHECK
    then
	echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s
    fi
    
    for MODE in HIGH STD; do
        typewalk --tags GHOST FLAT $CAM $MODE --dir $FLATDIR/ -n -o flat.list
        reduce --drpkg ghostdr @flat.list

	if $CHECK
	then
	    echo 'You can now check the reduction at this step.'
	    read -p "Press any key to continue... " -n1 -s
	fi
	
        typewalk --tags GHOST ARC $CAM $MODE --dir $ARCDIR/ -n -o arc.list
        reduce --drpkg ghostdr @arc.list 

	if $CHECK
	then
	    echo 'You can now check the reduction at this step.'
	    read -p "Press any key to continue... " -n1 -s
	fi

	###### NOTE THAT THE SINGLE SEEING IS BEING USED HERE INSTEAD OF BOTH ########
	###### The pipeline combines too many things. Assumed is also that any
	# 1.0 seeing files have also been removed, including the slitv ones.
	while read object <&3; do
            echo Reducing $object
            reduce --drpkg ghostdr $OBJDIR/$object
	    
	    if $CHECK
	    then
		echo 'You can now check the reduction at this step.'
		read -p "Press any key to continue... " -n1 -s
	    fi
	    
        done 3< <(
            typewalk --tags GHOST $CAM $MODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*$BINNING.*\.(fits|FITS)" -n \
                | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" \
                | grep '^\s' \
                | awk '{print $1}'
        )
    done
done
