#!/bin/bash

# make it so CTRL-C kills *all* subprocesses (doesn't require the user to press CTRL-C multiple times)
trap 'kill -s KILL -- -$$ 2>/dev/null' EXIT
trap 'exit' INT QUIT TERM

#Script to reduce all GHOST data using the recipe system.
#This must be ran from within where the GHOST data are kept, with default variables indicating where the various types of files are

# Start by setting up all the locations for the files. Change this for each case.
COREDIR=$PWD

CHECK=true

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

# typewalk --tags GHOST SLITV BIAS --dir $BIASDIR/ -n -o bias.list
# reduce --drpkg ghostdr @bias.list

# if $CHECK
# then
#     echo 'You can now check the reduction at this step.'
#     read -p "Press any key to continue... " -n1 -s
# fi    

# typewalk --tags GHOST SLITV DARK --dir $DARKDIR/ -n -o dark.list
# reduce --drpkg ghostdr @dark.list

# if $CHECK
# then
#     echo 'You can now check the reduction at this step.'
#     read -p "Press any key to continue... " -n1 -s
# fi

for mode in high std; do
    CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

    # typewalk --tags GHOST SLITV FLAT $CAPMODE --dir $FLATDIR/ -n -o flat.list
    # reduce --drpkg ghostdr @flat.list

    # if $CHECK
    # then
    # 	echo 'You can now check the reduction at this step.'
    # 	read -p "Press any key to continue... " -n1 -s
    # fi
    
    # typewalk --tags GHOST SLITV ARC $CAPMODE --dir $ARCDIR/ -n -o arc.list
    # reduce --drpkg ghostdr @arc.list

    # if $CHECK
    # then
    # 	echo 'You can now check the reduction at this step.'
    # 	read -p "Press any key to continue... " -n1 -s
    # fi
    
    while read object <&3; do
        echo Reducing $object
        reduce --drpkg ghostdr $object
	if $CHECK
	then
	    echo 'You can now check the reduction at this step.'
	    read -p "Press any key to continue... " -n1 -s
	fi
    done 3< <(
        typewalk --tags GHOST SLITV IMAGE $CAPMODE --dir $OBJDIR/ --filemask 'obj.*$SEEING.*\.(fits|FITS)' \
            -n -o tmp$$.list >& /dev/null && cat tmp$$.list | grep -v '^#'; rm tmp$$.list
    )
done

echo 'Now the spectrograph data'

for cam in red blue; do
    echo "Doing $cam images now"
    CAPCAM=`echo $cam | tr '[:lower:]' '[:upper:]'`

    typewalk --tags GHOST BIAS $CAPCAM --dir $BIASDIR/ --filemask '.*1x1.*\.(fits|FITS)' -n -o bias.list
    reduce --drpkg ghostdr @bias.list
    
    if [ $BINNING != '1x1' ]
    then
	typewalk --tags GHOST BIAS $CAPCAM --dir $BIASDIR/ --filemask '.*'$BINNING'.*\.(fits|FITS)' -n -o bias.list
	reduce --drpkg ghostdr @bias.list
    fi

    if $CHECK
    then
	echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s
    fi
    
    typewalk --tags GHOST DARK $CAPCAM --dir $DARKDIR/ -n -o dark.list
    reduce --drpkg ghostdr @dark.list
    
    if $CHECK
    then
	echo 'You can now check the reduction at this step.'
	read -p "Press any key to continue... " -n1 -s
    fi
    
    for mode in high std; do
	CAPMODE=`echo $mode | tr '[:lower:]' '[:upper:]'`

        typewalk --tags GHOST FLAT $CAPCAM $CAPMODE --dir $FLATDIR/ -n -o flat.list
        reduce --drpkg ghostdr @flat.list

	if $CHECK
	then
	    echo 'You can now check the reduction at this step.'
	    read -p "Press any key to continue... " -n1 -s
	fi
	
        typewalk --tags GHOST ARC $CAPCAM $CAPMODE --dir $ARCDIR/ -n -o arc.list
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
            reduce --drpkg ghostdr $object
	    
	    if $CHECK
	    then
		echo 'You can now check the reduction at this step.'
		read -p "Press any key to continue... " -n1 -s
	    fi
	    
        done 3< <(
            typewalk --tags GHOST $CAPCAM $CAPMODE --dir $OBJDIR/ --filemask "obj.*$SEEING.*$BINNING.*\.(fits|FITS)" \
                     -n -o tmp$$.list >& /dev/null && cat tmp$$.list | grep -v '^#'; rm tmp$$.list
        )
    done
done
