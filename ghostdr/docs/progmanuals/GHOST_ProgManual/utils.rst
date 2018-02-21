Utility Scripts
===============

This chapter describes various utility scripts provided with the Ghost DRS.
Some are for trouble-shooting while others are intended for use by Gemini
Observatory staff and the development team during the commissioning process.

``slit_viewer.py``
---------------------------------

The ``slit_viewer.py`` script takes a single file argument, an optional
``-save`` flag, and an optional ``-no_simult`` flag on the command line:

  usage: slit_viewer.py [-h] [-save] [-no_simult] path

The script clips the 2d profiles out of every slit viewer frame within the
file and displays them alongside 1d versions of the same.  The profiles are
shown in separate windows for each arm.  The profiles are arranged vertically,
with the 2d images amalgamated together and shown to the left of the 1d graphs:

   .. only:: latex

       .. image:: images/std-blue-with-cosmic.png
         :scale: 70
       .. image:: images/high-blue-with-arc-removed.png
         :scale: 70

   .. only:: html

       .. image:: images/std-blue-with-cosmic.png
         :scale: 65
       .. image:: images/high-blue-with-arc-removed.png
         :scale: 65

The displayed images and graphs allow one to easily see where, if any, cosmic
rays may have fallen on the fibres, and the relative levels of the background
between slit viewer frames.

The ``-save`` option saves the composite/amalgamated images to FITS files
adjacent to the source file (labelled ``_red.FITS`` and ``_blue.FITS``
respectively), while the ``-no_simult`` option zeros the simultaneous arc fibre
pixels in the profiles (both 2d and 1d) before generating/displaying the plots.
The ``-no_simult`` option is only valid for (only operates on) High resolution
OBJECT frames; otherwise it is ignored.


``input_locations.py``
----------------------

This is a utility script to remove complexity for the first part of the utility scripts.
It essentially contains variables with the locations for each file needed by the
various scripts and should be up-to-date in order for them to run successfully.

Given that the utility scripts are not meant to run pipeline primitives but instead
test parts of the pipeline by importing the core code itself, files must be given
as inputs. They, therefore, now require that this particular script is populated
correctly. It takes a mode, camera and user as inputs to a Files class, and then
picks up the default model files from the lookups directory and, depending on each user,
the base directory containing the raw files and calibrated files as well as inputs.

Please see the script itself for more information on what is needed. 


``ghost_slider_adjust.py``
---------------------------------

This script is a handy tool designed to allow users to visualise the xmod and
wavelength scale model superimposed with the data it is being fitted to.
It also allows the users to manually change the model to fit any new data and save
the result into a file.

This script allows users to select either the X model or the W model depending on
what they desire, and will superimpose a given initial model with a provided image,
typically convolution of the flat field with the slit profile (for the xmod) or an
arc frame (for the Wavelength model).

Users may define what files to use within the script itself and users are encouraged
to use the input_locations module for this purpose.

The scripts starts out gathering all necessary files, images and modules (for each
model case) and initiates the ghost class.

In any case, the core of this script is to run the ``ghost.manual_model_adjust`` function
which does all the work for this utility. For more information on this function see the
documentation for the polyspect module. 

In the case of the X model, the script will run the flat field convolution with the
flat image and then call the adjust function. Only the xparams is required here.

The W model will just call the function. In both cases the used can save the parameters to
a file. 

The 'percentage_variation' variable refers to the the amount that the sliders is allowed
to vary. For initial manipulation you may want to set this high and then very low for
finer adjustments.



``slit_rotation_fit.py``
----------------------------------

This script is designed to be used by the commissioning team to determine the
slit rotation model for the spectrograph once the first arc frame is taken.

All the file inputs are determined within the script and using the input_locations,
much like the ghost_slider_adjust.py
tool and the purpose is to determine the rotation of the slit as a function of order and
pixel along the orders.

This will produce a model similar to xmod and wmod which is used within the spectral format
function.

This tool will extract in 1D the arc frame and splits each order into sections (8 by default).
It then loops through every order and section taking a cross correlation between two extracted
profiles (objects 1 and 2 in std resolution mode; the object and sky in high resolution mode).
The trigonometric solution of this cross correlation and the pixel separation between the centroid
of each profile gives the angle for that section.

Angles outside of +/- 10 degrees or those sections with no flux are given 0 weighting to avoid
problems.

The resulting measured angles for all sections are then fitted in a least squares approach
to determine the best model for the slit rotation across the ccd.


``spatial_scale_fit.py``
----------------------------------

This script is designed to be used by the commissioning team to determine the
spatial scale model for the spectrograph once the first flat frame is taken.

All the file inputs are determined within the script, much like the ghost_slider_adjust.py
tool and the purpose is to determine the scale in the spatial direction of
the slit as a function of order and pixel along the orders.

This will produce a model similar to xmod and wmod which is used within the spectral format
function.

This tool takes a flat field frame and splits each order into sections (8 by default).
It then uses a range of pre-determined suitable scales (in slit microns per detector pixel) and
tests all scenarios by convolving the slit profile scaled (in a fixed manner) to each value
with the flat field using the pre-existing ``slit_flat_convolve`` function of the ghost class.

It then loops through every order and section taking the maximum of the convolution and corresponding
scale for each section, in order to determine what is the best scale for each case.

The resulting scales measured for all sections are then fitted in a least squares approach
to determine the best model for the spectrograph. This fit is weighted by the maximum convolution value
to ensure the edges of the chip where the flat flux is low are down weighted.


``quick_wl.py``
---------------

This script is here so that users can do a quick wavelength fitting test and inspect the results.

Files and models are imported in the standard way, and then the extracted arc flux is fed into the
``find_lines`` function, designed to look for the actual arc line positions. This function now contains
two inspection methods with boolean inputs: 'inspect', which will show the full image with blue crosses
over the initial position guess and red crosses over the found arc lines; 'plots' , if True, will show a
two panel plot with the extracted flux superimposed on the gaussian fit on the left, and a thumbnail of
the arc image on the right for each line in the finding routine. Going through these can take a lot of time
but it may be useful to determine what exactly is going wrong.

The 'before' fit and 'after' fit results are shown using the ``manual_model_adjust`` function.


``diagnostics.py``
------------------

This script is designed to be ran inside the data directory *after* the data has finished reducing.

This script can take optionally 1 or 2 command line inputs indicating mode and/or camera in any order.

E.g: ``diagnostics.py red high``
or:  ``diagnostics.py blue``

If any inputs are given, the script will inspect only those specified. Otherwise all 4 combinations
will be picked up and shown.

This will then show the X model superimposed on the convlution and raw flat field for inspection,
followed by a three panel plot containing the extracted arc flux with the arc line spectrum
superimposed for all 3 extracted objects.

The arc orders alternate in color between red and blue for convenience. The ideal arc spectrum is
plotted in green.

This allows the user to inspect the results of the reduction and ensure the pipeline is working
within the expected parameters.

