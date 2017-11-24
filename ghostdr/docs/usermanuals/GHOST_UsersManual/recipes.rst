.. recipes:

.. _GHOST_Recipes_and_Flows:

*****************
Recipes for GHOST
*****************

Setting up the development DRAGONS system
=========================================

For this discussion, we assume that you've already installed
a compatible Anaconda Python distribution, and done the setup necessary to
create a Gemini data reduction environment. For the purposes of this manual,
we assume that your environment is named ``geminidev``. Therefore, before doing
anything, you'll need to activate the environment::

    source activate geminidev

Furthermore, for the commands given below to work properly, you must:

 #. Add the following line to the top of ``DRAGONS/gempy/scripts/typewalk.py``::

        import ghost_instruments

 #. Add the following paths to your system ``PATH`` variable::

        DRAGONS/recipe_system/scripts
        GHOSTDR/scripts
        DRAGONS/gempy/scripts

 #. Add the following paths to your environment ``PYTHONPATH`` variable::

        DRAGONS
        GHOSTDR/simulator/pyghost
        GHOSTDR

Typical Processing Flows
========================

This section describes how to process simulated GHOST data. For information on
how to use the simulator to generate simulated GHOST data, please refer to
its documentation: :any:`ghostsim`.

Setting up the calibrations system
----------------------------------

At present, we are using a beta version of the Gemini local calibration
manager. Assuming that you have access to the current version of this code
(at the time of writing, ``GeminiCalMgr-0.9.9.6-ghost``), you need to take
the following steps to prepare the calibration manager for use:

#. Activate your ``geminidev`` Anaconda environment;
#. Install the calibration manager::

    python /path/to/GeminiCalMgr-0.9.9.6-ghost/setup.py install

#. Initialize your local database, something like this (on a Unix system; for
   MacOSX or Windows, use an appropriate file path for your calibration file)::

    mkdir /home/you/calmgrdb
    mkdir /home/you/.geminidr
    cat -> /home/you/.geminidr/rsys.cfg <<-HERE
        [calibs]
        standalone = True
        database_dir = /home/you/calmgrdb/cal_manager.db
    HERE
    reduce_db.py init -v -w

#. Apply the following GHOST-related patches to the calibration system code:

   - ``/path/to/GeminiCalMgr-0.9.9.6-ghost/src/cal/calibration_ghost``:

        - Remove ``Ghost.nodandshuffle`` from around line 183;
        - Add the following at around line 345::

            def processed_slitflat(self, howmany=None):
                return self.flat(True, howmany)


   - ``/path/to/GeminiCalMgr-0.9.9.6-ghost/src/orm/calibration_ghost``:
        - Add the following around line 15::

            RESOLUTIONS = ['std', 'high']
            RESOLUTION_ENUM = Enum(*RESOLUTIONS, name='ghost_resolution')

        - Add the following around line 35::

            res_mode = Column(RESOLUTION_ENUM, index=True)

        - Add the following around line 67::

            resolution = ad.res_mode()
            if resolution in RESOLUTIONS:
                self.res_mode = resolution

#. Deploy the changes you just made to the calibration system::

    python /path/to/GeminiCalMgr-0.9.9.6-ghost/setup.py install

DRAGONS does not currently automatically send its output calibration files to
the GeminiCalMgr. You will have to do this manually after each step, e.g.::

    reduce_db.py add calibrations/processed_thing/my_processed_thing.fits

where ``thing`` is ``bias``, ``flat``, ``dark``, etc.

.. note::
    The calibration system knows about which context (quality assurance, ``qa``,
    or science quality, ``sq``) a calibrator was generated in, and will only
    send you back that calibrator if you are in the same context (so, you
    cannot retrieve ``qa`` calibrators if you are attempting an ``sq``
    reduction, or vice versa). The default context is ``qa``.

    You can switch contexts using the ``--context`` flag to ``reduce``.

Bulk-reducing simulated data
----------------------------

To simplify the end-to-end testing of :any:`ghostdr`, Joao Bento has
contributed a bash script to run the entire reduction sequence in one command::

    ./GHOSTDR/utils/reduce_all.sh

The script will run reductions for all combinations of spectrograph arm and
resolution in sequence. The script will pause between each distinct
reduction step (e.g. bias reduction, dark reduction, etc.) to allow you to
review the output. Alternatively, you could choose to run through the steps
manually yourself, as described below.

Data Reduction Steps
--------------------

The reduction of each component of a GHOSTDR observing package (bias, dark,
flat, etc.) can be broken down into three parts:

``typewalk``
++++++++++++

.. note::
    For reducing a single file, you don't need to use the ``typewalk``
    utility.

The ``typewalk`` utility is used for generating lists of files to reduce
together. This may be because the list of files will eventually require
stacking, or simply as a convenience for reducing a number of data frames
with a single command.

The most common usage for ``typewalk`` is to generate a list of files
with matching :module:`AstroData` tags. For example, to generate a list
of all files in the current directory which are red camera biases with 2x4
binning, and write this list out to a text file called
``bias.1x1.red.list``, use the following::

    typewalk --tags GHOST BIAS RED 1x1 -o bias.1x1.red.list

There are several other options available (e.g. using a regex filemask to
further restrict the files you're considering) -- type ``typewalk --help`` to
see these options.

``reduce``
++++++++++

The ``reduce`` command is part of :any:`DRAGONS`, and works
in a similar fashion to an old ``IRAF`` call. Please see the :any:`DRAGONS`
documentation for more detail. However, there are two important options to
take note of for development GHOST reduction::

    reduce --drpkg ghostdr @bias.1x1.red.list

The option ``--drpkg ghostdr`` tells ``reduce`` to import the ``ghost``
data reduction package, in addition to the standard :any:`DRAGONS`. This
will not be required in production, and ``ghostdr`` will be incorporated
into :any:`DRAGONS` by Gemini.

The ``@`` modifier tells ``reduce`` that the input file is in fact a list,
and should be broken apart for reduction. If you were only passing a single
FITS file to ``reduce``, you would leave the ``@`` modifier off.

``reduce_db.py``
++++++++++++++++

The current iteration of the local calibration manager has no ability to
automatically detect when a new calibrator has appeared in the
``calibrations/`` directory. Therefore, you will need to manually load
your calibrators into the system::

    reduce_db.py add calibrations/processed_bias/bias_1_1x1_red_bias.fits

The ``reduce_db.py remove`` command has the same syntax, and can be used to
remove files from the database. This is useful if your original calibrator
has been superseded, or you've accidentally added a file to the database you
shouldn't have (e.g. a rebinned dark or flat). To see all the files
currently referenced in the database, use::

    reduce_db.py list

.. _reducing-slit-viewing-images:

Reducing Slit Viewing Images
----------------------------

The first step in reduction is to create slit viewer frames
(which, when applied, remove cosmic rays and
compute the mean exposure epoch).  The first step, computing the slit bias
calibrator, may be skipped in favour of simply pointing to a single slit bias
frame
(with tags ``SLITV`` and ``BIAS``).  Or, follow these steps to produce one by
stacking multiple frames together::

    typewalk --tags GHOST BIAS SLITV --dir <path_to>/data_folder -o slit.bias.list
    reduce --drpkg ghostdr @slit.bias.list
    reduce_db.py add calibrations/processed_bias/your_red_SLIT_bias.fits

The next step is to generate the dark calibrator.  Follow these steps to produce
one::

    typewalk --tags GHOST SLITV DARK --dir <path_to>/data_folder -o slit.dark.list
    reduce --drpkg ghostdr @slit.dark.list
    reduce_db.py add calibrations/processed_dark/your_red_SLIT_dark.fits

Now generate the flat calibrator.  For this you will now need to specify an
additional type to ``typewalk`` that identifies the resolution of the data that
you wish to process (as mixing resolutions would be nonsensical).  Follow these
steps as an example::

    typewalk --tags GHOST SLITV FLAT STD --dir <path_to>/data_folder -o slit.flat.std.list
    reduce --drpkg ghostdr @slit.flat.std.list
    reduce_db.py add calibrations/processed_slitflat/your_red_SLIT_slitflat.fits

Though not (yet) used in our final object reduction, you can also produce a
master arc frame::

    typewalk --tags GHOST SLITV ARC HIGH --dir <path_to>/data_folder -o slit.arc.std.list
    reduce --drpkg ghostdr @slit.arc.std.list
    reduce_db.py add calibrations/processed_arc/your_red_SLIT_arc.fits

The final step is to use all of the above calibrators (except the arc) in a call
to ``reduce`` a set of slit viewer images taken concurrently with a science
frame, usually found in files named like ``obj95_1.0_high_SLIT.fits`` (following
this convention: ``obj{exptime}_{seeing}_{resolution}_SLIT.fits``).
This informs the reduction framework to run the
``makeProcessedSlit`` GHOST recipe on them.  Run the reduction as follows::

    reduce --drpkg ghostdr <path_to>/data_folder/obj95_1.0_high_SLIT.fits

Generating a Bias Calibration frame
-----------------------------------

Once you have a few biases of the same arm to work with, generate a file list
using the ``typewalk`` utility.  The following command assumes you have
generated several red arm biases with a 1x1 binning::

    typewalk --tags GHOST BIAS RED 1x1 --dir <path_to>/data_folder -o bias.1x1.red.list

The ``--dir`` argument can be omitted if you are already within the folder
containing the data.

.. warning::
    Make sure you've made the necessary changes to the ``typewalk.py`` script!

Now you are ready to generate a bias calibration frame.  The following command
(which runs the ``makeProcessedBiasG`` Gemini recipe behind the scenes) will
stack the bias frames in listed ``bias_red.list`` and store the finished bias
calibration in ``calibrations/processed_bias/``::

    reduce --drpkg ghostdr @<path_to>/bias.1x1.red.list
    reduce_db.py add calibrations/processed_bias/your_red_bias.fits

Don't forget the @ character in this line, e.g. if <path_to> is ``data`` then
this command should be ``reduce @data/bias.list``. The @ parameter is a legacy
from IRAF, and tells ``reduce`` that you're passing a list of filenames instead
of a data file.

.. note::
    This example uses 1x1 binned data. If you are reducing data in another
    binning mode, you will need to reduce the biases of that binning mode,
    *as well as* the standard 1x1 binned biases. This is because darks, arcs and
    flats are always taken at 1x1 binning, so require reduced 1x1 binned
    biases to be reduced correctly.


The ``--drpkg ghostdr`` flag tells the recipe system it should attempt to import
from the ``ghostdr`` folder (which should now be on your ``PYTHONPATH``), in
addition to the standard ``DRAGONS`` system. In production, this flag will
be unnecessary; ``ghostdr`` will come as part of ``DRAGONS``.

This code call will place a file named something like ``bias_1_red_bias.fits``
in the
``calibrations/processed_bias`` directory of your present working directory.
This file will then be added to the calibrations directory by the
``reduce_db.py`` script call.

.. note::
    The final name of stacked frames (of which your bias is one) depends on
    which input file was queued up to be stacked first. This, in turn,
    depends on the output of an :any:`os.listdir` call, which returns files
    in disk order, *not* name order (like the ``ls`` system command does on
    Unix). Therefore, it cannot be guaranteed that your stacked bias file name
    will be ``bias_1_red_bias.fits`` - among other things, the number in the
    middle may be different.

The whole process behind Gemini's ``makeProcessedBias`` recipe is documented in
the following flowchart (thanks Kathleen Labrie):

.. only:: latex

    .. image:: images/biasCalibration.png
      :scale: 70

.. only:: html

    .. image:: images/biasCalibration.png
      :scale: 45

Generating a Dark Calibration Frame
-----------------------------------

The procedure for generating a dark calibration frame is broadly similar to
making a bias calibration frame. However, the tags to be passed to ``typewalk``
should be ``DARK`` instead of ``GHOST_BIAS`` (in addition to the
necessary ``RED``/``BLUE`` tag)::

    typewalk --tags DARK RED --dir <path_to>/data_folder -o dark.red.list

The dark frames may then be reduced by invoking::

    reduce --drpkg ghostdr @<path_to>/dark.red.list
    reduce_db.py add calibrations/processed_dark/your_red_dark.fits

Make sure you've run ``prepare_data.py`` over your data directory before
attempting this step, otherwise the reduction system will not be able to locate
your previous bias calibration.

The whole process behind Gemini's ``makeProcessedDark`` recipe is documented in
the following flowchart (thanks Kathleen Labrie):

.. only:: latex

  .. image:: images/darkCalibration.png
    :scale: 70

.. only:: html

  .. image:: images/darkCalibration.png
    :scale: 45


Generating a Flat Calibration Frame
-----------------------------------

.. warning:: You *must* have performed a full slit viewer reduction before
             attempting to make a flat calibrator. See
             :ref:`reducing-slit-viewing-images` for details.

The procedure for generating a flat field calibration frame is similar to
creating a dark or bias, although you have to ``typewalk`` over FLAT files
instead. You also need to specify an instrument resolution for the first time,
e.g.::

    typewalk --types FLAT GHOST STD RED --dir <path_to>/data_folder -o flat.red.std.list

A simple call to ``reduce`` once again processes the list of flats::

    reduce --drpkg ghostdr @<path_to>/flat.red.std.list
    reduce_db.py add calibrations/processed_flat/your_red_flat.fits

After the flat field has been created, the spectrograph apertures are fit using
a ``polyfit`` approach. ``DRAGONS`` will read in the appropriate aperture
model from the ``lookups`` system, fit it to the flat field, and store the
resulting model in the calibrations system.

The selection of the appropriate ``polyfit`` model to start with is
determined by the spectrograph arm, resolution, and the date the observations
are made on. Ideally, there will only be one model per arm and resolution
combination; however, spectrograph maintenance (i.e. dis- and re-assembly) may
result in the model changing at a specific point in time. Therefore, the
RecipeSystem *should* (see below) automatically choose the most recent
applicable model for the dataset being considered.

.. note:: Date-based model selection is currently not implemented - instead,
          only a single model is provided for each arm/resolution combination.
          This is sufficient for testing involving the simulator data.
          Date-based selection will be implemented soon.

The process behind ``makeProcessedFlatG`` is summarized in the following
flowchart (thanks Kathleen Labrie):

.. only:: latex

    .. image:: images/flatCalibration.png
      :scale: 70

.. only:: html

    .. image:: images/flatCalibration.png
      :scale: 45

.. note:: This is the originally-envisaged implementation of
          ``makeProcessedFlatG``. It has since been decided that Gemini will
          guarantee that Gemini Observatory will always take at least three
          flat fields per arm per observation, which means that
          ``rejectCosmicRays`` is not required; ``stackFrames`` will remove
          almost all cosmic rays.


Generating an Arc Calibration Frame
-----------------------------------

.. warning:: You *must* have performed a full slit viewer reduction before
             attempting to make an arc calibrator - the results of the slit
             flat and slit image reduction are required to make the profile
             extraction and subsequent wavelength fitting work. See
             :ref:`reducing-slit-viewing-images` for details.

Making an arc calibration frame is similar to the previous calibration steps.
The correct tag to ``typewalk`` across is ``ARC``::

    typewalk --tags GHOST ARC RED STD --dir <path_to>/data_folder -o arc.red.std.list

Then, the following command reduces the arcs::

    reduce --drpkg ghostdr @<path_to>/arc.red.std.list
    reduce_db.py add calibrations/processed_arc/your_red_arc.fits

Arc reduction not only generates a reduced arc image and places it in the
calibrations directory, but also uses the ``polyfit`` module to extract the
flux profiles of the object/sky fibres in the input image. It then uses this
fit, and a line set stored in the RecipeSystem lookups system, to make a
wavelength fit to the arc image. This fit is also stored in the calibrations
directory/system.


Reducing an Object frame (Spectra)
----------------------------------

The GHOST simulator produces object spectra frames like
``obj95_1.0_high_red.fits`` whose names follow this convention:
``obj{exptime}_{seeing}_{resolution}_{arm}.fits``. If you run ``typewalk`` on
the folder containing these, you'll see that they are identified as having the
tag ``SPECT``, but none of the further tags we've encountered already (e.g.
``BIAS``, ``DARK``, etc.)::

    typewalk --dir <path_to>/data_folder

This informs the reduction framework to run the ``reduceG`` GHOST recipe on
them. which should run to at least the ``flatCorrect`` step now that you
have dark and bias calibration frames (for the moment, we have commented the
remaining steps out of the ``reduceG`` recipe so it will complete
successfully)::

    reduce --drpkg ghostdr <path_to>/data_folder/obj95_1.0_std_1x1_red.fits

This produces a ``obj95_1.0_std_1x1_red_flatCorrected.fits`` (or similar) file,
a bias, dark and flat corrected GHOST spectrum frame.

.. warning:: The primitive ``rejectCosmicRays`` would normally be called as
             part of ``reduceG``, after the ``darkCorrect`` step. It is
             currently commented out - the underlying LACosmic algorithm is
             working, but aperture removal/re-instatement is required to avoid
             accidentally flagging spectral peaks and the edges of orders as
             cosmic rays, and this has yet to be implemented.


Other Processing Flows
======================
include scientific flow charts, include associated recipes
