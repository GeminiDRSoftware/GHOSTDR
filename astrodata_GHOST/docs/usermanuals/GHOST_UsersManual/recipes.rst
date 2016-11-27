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


Reducing an Object frame (Spectra)
----------------------------------

The GHOST simulator produces object spectra frames like
``obj100_1.0_std_red.fits`` whose names follow this convention:
``obj{exptime}_{seeing}_{resolution}_{arm}.fits``. If you run ``typewalk`` on
the folder containing these, you'll see that they are identified as
``GHOST_SPECT``::

    typewalk --dir <path_to>/data_folder

This informs the reduction framework to run the ``reduceG`` GHOST recipe on
them, which should run to at least the ``darkCorrection`` step now that you
have dark and bias calibration frames (for the moment, we have commented the
remaining steps out of the ``reduceG`` recipe so it will complete
successfully::

    reduce <path_to>/data_folder/obj100_1.0_std_red.fits

The above command will fail due to the faulty calibrations lookup. Again, we
need to use the ``--override_cal`` option::

    reduce <path_to>/data_folder/obj100_red.fits --override_cal
    processed_bias:calibrations/storedcals/bias_0_red_bias.fits
    processed_dark:calibrations/storedcals/dark100_0_red_dark.fits

This produces a ``obj100_1.0_std_red_darkCorrected.fits`` (or similar) file, a
bias and dark corrected GHOST spectrum frame.


Finding apertures
-----------------

This is the first iteration of this documentation. At this stage this process
is done purely within python. 

The principle behind this process is that the format of the spectral orders 
and the wavelength scale can be modelled uniquely using polynomials of
polynomials. e.g. A series of polynomials as a function of order number are 
combined in a polynomial fashion as a function of y position on the chip. 

The polyfit module is used at this stage and requires knowledge of the 
spectrograph arm, mode and a reduced flat field image (tested only with 
flats from the simulator). 

Usage follows::

  import polyfit
  import astropy.io.fits as pyfits
  ghost_format = polyfit.ghost.Arm('red',mode='high')

At this stage it is important to have some existing file or initial guess
array for the polynomial model. By default, these are found in a 
``../data/ghost/<arm>/<mode>/`` folder alongisde the polyfit module but this 
can be changed using the ``ghost_format.model_location`` variable. 

After acquiring the flat field data::

  flat_file = "location_to_flat/flatfield_frame.fits"
  flat_data = pyfits.getdata(flat_file)

a convolution map is required. This is done so that, irrespective of the 
number of fibers per order, the model is adjusted with respect to an 
equivalent map that has maxima where the middle of the order lies. 

Either a supplied model of the slit profile (from the slit viewer) or a 
default uniform illumination profile is convolved with every column of the 
flat field image along the spatial direction, resulting in a series of 
images that match the centers of the orders. 

This is then fed into the ``adjust_model`` function for visual inspection
of the initial model::

  flat_conv=ghost_format.slit_flat_convolve(flat_data)
  ghost_format.adjust_model(flat_conv,convolve=False,percentage_variation=10)

The ``percentage_variation`` refers to the percentage range of values that each
parameter is allowed to be varied by the matplotlib slider widgets. 

The ``Submit`` button saves the current version of the model onto the default 
location. 

Once the model is close, the model can be fitted::

  ghost_format.fit_x_to_image(flat_conv,decrease_dim=8,inspect=True)

This function takes the convolution map and adjusts the model to the local
maximum along each order. The ``inspect`` parameter set to True displays the
result of the fit and, if it is satisfactory, pressing the ``submit`` button
again will save it. 

At this point, the input to the flux extraction code is available for testing.


Other Processing Flows
======================
include scientific flow charts, include associated recipes

