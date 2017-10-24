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

 #. There is currently an issue in the :any:`stackFrames` primitive which causes
    unexpected behaviour when using the cheat calibration setup. To rectify
    this, add the following below line 58 of
    ``DRAGONS/geminidr/core/primitives_bookkeeping.py``::

        for i in self.stacks: self.stacks[i].sort()

    This issue will be rectified once the full calibration manager is up and
    running.
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

At present, we do not have a functioning local calibration system. In the
meantime, Chris Simpson from Gemini has provided us with a script that creates
an association dictionary between the simulator observation ID, file name and
the name of the expected calibrator(s). It should be run in the same directory
as your data files::

    python prepare_data.py

This script will generate a dictionary keyed by observation ID and required
calibration type, the value of which is the name of the expected calibration
file (which will eventually be written into the appropriate
``calibrations/processed_*`` directory.

The output of this script is pickled into the file ``.reducecache/calindex.pkl``
within your working directory. You can inspect it by loading it into a
(i)Python shell::

    import pickle
    with open('./reducecache/calindex.pkl', 'r') as fileobj:
        caldict = pickle.load(fileobj)

.. note::
    The expected calibrator names are hard-coded into ``prepare_data.py``,
    and assume that you haven't changed the data file names from what the
    simulator generated. If you do need/want to change the file names for
    whatever reason, you will need to modify ``prepare_data.py``, or load
    the calibrations dictionary, edit it yourself, and re-pickle it to disk.

Bulk-reducing simulated data
----------------------------

To simplify the end-to-end testing of :any:``ghostdr``, Joao Bento has
contributed a bash script to run the entire reduction sequence in one command::

    ./GHOSTDR/utils/reduce_all.sh

The script will run reductions for all combinations of spectrograph arm and
resolution in sequence. The script will pause between each distinct
reduction step (e.g. bias reduction, dark reduction, etc.) to allow you to
review the output. Alternatively, you could choose to run through the steps
manually yourself, as described below.

Generating a Bias Calibration frame
-----------------------------------

Once you have a few biases of the same arm to work with, generate a file list
using the ``typewalk`` utility.  The following command assumes you have
generated several red arm biases::

    typewalk --tags BIAS RED --dir <path_to>/data_folder -o bias_red.list

The ``--dir`` argument can be omitted if you are already within the folder
containing the data.

.. warning::
    Make sure you've made the necessary changes to the ``typewalk.py`` script!

Now you are ready to generate a bias calibration frame.  The following command
(which runs the ``makeProcessedBiasG`` Gemini recipe behind the scenes) will
stack the bias frames in listed ``bias_red.list`` and store the finished bias
calibration in ``calibrations/processed_bias/``::

    reduce --drpkg ghostdr @<path_to>/bias.list

Don't forget the @ character in this line, e.g. if <path_to> is ``data`` then
this command should be ``reduce @data/bias.list``. The @ parameter is a legacy
from IRAF, and tells ``reduce`` that you're passing a list of filenames instead
of a data file.

The ``--drpkg ghostdr`` flag tells the recipe system it should attempt to import
from the ``ghostdr`` folder (which should now be on your ``PYTHONPATH``), in
addition to the standard ``DRAGONS`` system. In production, this flag will
be unnecessary; ``ghostdr`` will come as part of ``DRAGONS``.

This code call will place a file named ``bias_1_red_bias.fits`` in the
``calibrations/processed_bias`` directory of your present working directory.

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

    typewalk --tags DARK RED --dir <path_to>/data_folder -o dark_red.list

The dark frames may then be reduced by invoking::

    reduce --drpkg ghostdr @<path_to>/dark.list

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

    typewalk --types FLAT GHOST HIGH --dir <path_to>/data_folder -o flat_red.list

A simple call to ``reduce`` once again processes the list of flats::

    reduce --drpkg ghostdr @<path_to>/flat.list

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

    typewalk --tags ARC RED HIGH --dir <path_to>/data_folder -o arc_red.list

Then, the following command reduces the arcs:

    reduce --drpkg ghostdr @<path_to>/arc.list

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

    reduce --drpkg ghostdr <path_to>/data_folder/obj95_1.0_high_red.fits

This produces a ``obj95_1.0_high_1x1_red_flatCorrected.fits`` (or similar) file, a
bias, dark and flat corrected GHOST spectrum frame.

.. warning:: The primitive ``rejectCosmicRays`` would normally be called as
             part of ``reduceG``, after the ``darkCorrect`` step. It is
             currently commented out - the underlying LACosmic algorithm is
             working, but aperture removal/re-instatement is required to avoid
             accidentally flagging spectral peaks and the edges of orders as
             cosmic rays, and this has yet to be implemented.

.. _reducing-slit-viewing-images:

Reducing Slit Viewing Images
----------------------------

Reducing slit viewer images is very similar to reducing standard images,
including steps to generate bias, dark and flat calibration frames, plus a
final step to process the slit viewer frames (which removes cosmic rays and
computes the mean exposure epoch).  The first step, computing the bias
calibrator, may be skipped in favour of simply pointing to a slit bias frame
(with tags ``SLITV`` and ``BIAS``).  Or, follow these steps to produce one by
stacking multiple frames together::

    typewalk --tags BIAS SLITV --dir <path_to>/data_folder -o slit_bias.list
    reduce --drpkg ghostdr @slit_bias.list

The next step is to generate the dark calibrator.  Follow these steps to produce
one::

    typewalk --tags SLITV DARK --dir <path_to>/data_folder -o slit_dark.list
    reduce --drpkg ghostdr @slit_dark.list

Now generate the flat calibrator.  For this you will now need to specify an
additional type to ``typewalk`` that identifies the resolution of the data that
you wish to process (as mixing resolutions would be nonsensical).  Follow these
steps as an example::

    typewalk --tags SLITV FLAT HIGH --dir <path_to>/data_folder -o slit_flat_high.list
    reduce --drpkg ghostdr @slit_flat_high.list

Though not (yet) used in our final object reduction, you can also produce a
master arc frame::

    typewalk --tags SLITV ARC HIGH --dir <path_to>/data_folder -o slit_arc_high.list
    reduce --drpkg ghostdr @slit_arc_high.list

The final step is to use all of the above calibrators (except the arc) in a call
to ``reduce`` a set of slit viewer images taken concurrently with a science
frame, usually found in files named like ``obj95_1.0_high_SLIT.fits`` (following
this convention: ``obj{exptime}_{seeing}_{resolution}_SLIT.fits``).
This informs the reduction framework to run the
``makeProcessedSlit`` GHOST recipe on them.  Run the reduction as follows::

    reduce --drpkg ghostdr <path_to>/data_folder/obj95_1.0_high_SLIT.fits


Other Processing Flows
======================
include scientific flow charts, include associated recipes
