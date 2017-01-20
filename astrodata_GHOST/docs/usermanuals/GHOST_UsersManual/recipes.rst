.. recipes:

.. _GHOST_Recipes_and_Flows:

*****************
Recipes for GHOST
*****************

Typical Processing Flows
========================

Here we review some of the typical processing workflows for GHOST data
reduction. For this discussion, it is assumed you've already installed the
latest Ureka package and made a local clone of the ``ghostdr`` Hg repository.
(For illustrative purposes, the text below assumes the clone's working copy
root is ``wc/``.)

Furthermore, for the commands given below to work properly, you must:

 #. initialize the Ureka environment: ``ur_setup``
 #. create a symlink named ``wc/externals/gemini_python/astrodata_GHOST``
    pointing to ``wc/astrodata_GHOST``,
 #. add ``wc/externals/gemini_python`` to the beginning of your ``PYTHONPATH``,
    and
 #. add ``wc/externals/gemini_python/astrodata/scripts`` and
    ``wc/externals/
    gemini_python/recipe_system/apps`` to the beginning of your
    ``PATH``

Generating a Bias Calibration frame
-----------------------------------

To generate a bias calibration frame you need 2 or more GHOST bias frames from
the same arm.  Until the instrument is live, you can use the GHOST simulator to
generate this data.  Its ``testsim.py`` script will create several types of
frames, including 3 bias frames for each arm, 3 darks for each arm, and 3 flats
for each arm and resolution combination. You can comment out the generation of
non-bias frame types to speed things up.

Once you have a few biases of the same arm to work with, generate a file list
using the ``typewalk`` utility.  The following command assumes you have
generated several red arm biases (if you don't specify either ``GHOST_RED`` or
``GHOST_BLUE``, you may get mixed red and blue frames which don't stack well!)::

    typewalk --types GHOST_BIAS GHOST_RED --dir <path_to>/data_folder -o
    bias.list

Now you are ready to generate a bias calibration frame.  The following command
(which runs the ``makeProcessedBiasG`` Gemini recipe behind the scenes) will
stack the bias frames in listed ``bias.list`` and store the finished bias
calibration in ``calibrations/storedcals/``::

    reduce @<path_to>/bias.list

Don't forget the @ character in this line, e.g. if <path_to> is ``data`` then
this command should be ``reduce @data/bias.list``. The @ parameter is a legacy
from IRAF, and tells ``reduce`` that you're passing a list of filenames instead
of a data file.
This code call will place a file named ``bias_0_red_bias.fits`` in the
``calibrations/storedcals`` directory of your present working directory.

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
making a bias calibration frame. However, the type to be passed to ``typewalk``
should be ``GHOST_DARK`` instead of ``GHOST_BIAS`` (in addition to the
necessary ``GHOST_RED``/``GHOST_BLUE`` type)::

    typewalk --types GHOST_DARK GHOST_RED --dir <path_to>/data_folder -o
    dark.list

Assuming ``typewalk`` has output your list of dark frames to ``dark.list``,
attempting to run::

    reduce @<path_to>/dark.list

will fail. This is because the framework cannot currently find calibrations
stored on disk (it uses a much more complicated lookup scheme).  The workaround
for the time being is to force it to look on disk in a particular area using the
``--override_cal`` option::

    reduce @<path_to>/dark.list  --override_cal
    processed_bias:calibrations/storedcals/bias_0_red_bias.fits

(Depending on your specific bias.list contents, your bias calibration under
your calibrations/storedcals directory may have a different name, so double-
check.) This command will place a file ``dark100_0_red_dark.fits`` into the
``calibrations/storedcals`` directory.

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

The procedure for generating a flat field calibration frame is similar to
creating a dark or bias, although you have to ``typewalk`` over GHOST_FLAT files
instead, e.g.::

    typewalk --types GHOST_FLAT GHOST_RED GHOST_HIGH --dir <path_to>/data_folder
    -o flat.list

(Note this is the first place where we have to explicitly specify the
resolution mode/type of the object file we ultimately intend to reduce.)
Then, when you call ``reduce`` on the ``flat.list``, you must provide both
the bias and flat file path explicitly::

    reduce @<path_to>/flat.list  --override_cal
    processed_bias:calibrations/storedcals/bias_0_red_bias.fits
    processed_dark:calibrations/storedcals/dark_0_red_dark.fits

(or whatever the filename of the processed dark turns out to be).

After the flat field has been created, the spectrograph apertures are fit using
a ``polyfit`` approach. The RecipeSystem will read in the appropriate aperture
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


Reducing an Object frame (Spectra)
----------------------------------

The GHOST simulator produces object spectra frames like
``obj100_1.0_std_red.fits`` whose names follow this convention:
``obj{exptime}_{seeing}_{resolution}_{arm}.fits``. If you run ``typewalk`` on
the folder containing these, you'll see that they are identified as
``GHOST_OBJECT``::

    typewalk --dir <path_to>/data_folder

This informs the reduction framework to run the ``reduceG`` GHOST recipe on
them. which should run to at least the ``flatCorrect`` step now that you
have dark and bias calibration frames (for the moment, we have commented the
remaining steps out of the ``reduceG`` recipe so it will complete
successfully)::

    reduce <path_to>/data_folder/obj100_1.0_high_red.fits

The above command will fail due to the faulty calibrations lookup. Again, we
need to use the ``--override_cal`` option::

    reduce <path_to>/data_folder/obj100_1.0_high_red.fits --override_cal
    processed_bias:calibrations/storedcals/bias_0_red_bias.fits
    processed_dark:calibrations/storedcals/dark100_0_red_dark.fits
    processed_flat:calibrations/storedcals/flat100_std_0_red_flat.fits

This produces a ``obj100_1.0_high_red_flatCorrected.fits`` (or similar) file, a
bias, dark and flat corrected GHOST spectrum frame.

.. warning:: The primitive ``rejectCosmicRays`` would normally be called as
             part of ``reduceG``, after the ``darkCorrect`` step. It is
             currently commented out - the underlying LACosmic algorithm is
             working, but aperture removal/re-instatement is required to avoid
             accidentally flagging spectral peaks and the edges of orders as
             cosmic rays, and this has yet to be implemented.


Reducing Slit Viewing Images
----------------------------

Reducing slit viewer images is very similar to reducing standard images,
including steps to generate bias and dark calibration frames, plus a final step
to process the slit viewer frames (which removes cosmic rays and computes the
mean exposure epoch).  The first step, computing the bias calibrator, may be
skipped in favour of simply pointing to a slit bias frame (of type
``GHOST_SLITV_BIAS``).  Or, follow these steps to produce one by stacking
multiple frames together::

    typewalk --types GHOST_SLITV_BIAS --dir <path_to>/data_folder -o
    slit_bias.list

    reduce @<path_to>/slit_bias.list

The next step is to generate the dark calibrator.  Follow these steps to produce
one::

    typewalk --types GHOST_SLITV_DARK --dir <path_to>/data_folder -o
    slit_dark.list

    reduce @<path_to>/slit_dark.list --override_cal
    processed_bias:calibrations/storedcals/bias_1_SLIT_bias.fits

A step that is not currently implemented (but will be for the next release) is
to generate the slit mask calibration as a byproduct of the normal flat image
processing.  When this is done, you will point to this calibrator using the
``--override_cal processed_slitmask:...`` option to the ``reduce`` script
invocation below.

The final step is to use all of the above calibrators in a call to ``reduce`` a
set of slit viewer images taken concurrently with a science frame, usually found
in files named as ``obj100_1.0_std_SLIT.fits`` (following this convention:
``obj{exptime}_{seeing}_{resolution}_SLIT.fits``). If you run ``typewalk`` on
the folder containing these, you'll see that they are identified as
``GHOST_SLITV_IMAGE``::

    typewalk --dir <path_to>/data_folder

This informs the reduction framework to run the ``makeProcessedSlitG`` GHOST
recipe on them.  Run the reduction as follows::

    reduce <path_to>/data_folder/obj100_1.0_std_SLIT.fits --override_cal
    processed_bias:calibrations/storedcals/bias_1_SLIT_bias.fits
    processed_dark:calibrations/storedcals/dark100_1_SLIT_dark.fits


Other Processing Flows
======================
include scientific flow charts, include associated recipes
