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



``ghost_slider_adjust.py``
---------------------------------

This script is a handy tool designed to allow users to visualise the xmod and
wavelength scale model superimposed with the data it is being fitted to.
It also allows the users to manually change the model to fit any new data and save
the result into a file.

This script allows users to select either the X model or the W model depending on
what they desire, and will superimpose a given initial model with either the
convolution of the flat field with the slit profile (for the xmod) or an arc frame
(for the Wavelength model).

Users may define what files to use within the script itself. Too many inputs would
be required and thus command line inputs are not available.

(EXAMPLE IMAGES TO FOLLOW)
..

The inner workings of this script will be described in a later release.


``slit_rotation_fit.py``
----------------------------------

This script is designed to be used by the commissioning team to determine the
slit rotation model for the spectrograph once the first arc frame is taken.

All the file inputs are determined within the script, much like the ghost_slider_adjust.py
tool and the purpose is to determine the rotation of the slit as a function of order and
pixel along the orders.

This will produce a model similar to xmod and wmod which is used within the spectral format
function.

This tool takes an 1D extracted arc frame and splits each order into sections (8 by default).
It then loops through every order and section taking a cross correlation between two extracted
profiles (objects 1 and 2 in std resolution mode; the object and sky in high resolution mode).
The trigonometric solution of this cross correlation and the pixel separation between the centroid
of each profile gives the angle for that section.
Angles outside of +/- 10 degrees or those sections with no flux are given 0 weighting.

The resulting measured angles for all sections are then fitted in a least squares approach
to determine the best model for the slit rotation across the ccd.



``spatial_scale_fit.py``
----------------------------------

This script is designed to be used by the commissioning team to determine the
spatial scale model for the spectrograph once the first flat frame is taken.

All the file inputs are determined within the script, much like the ghost_slider_adjust.py
tool and the purpose is to determine the rotation of scale in the spatial direction of
the slit as a function of order and pixel along the orders.

This will produce a model similar to xmod and wmod which is used within the spectral format
function.

This tool takes a flat field frame and splits each order into sections (8 by default).
It then uses a range of pre-determined suitable scales (in slit microns per detector pixel) and
tests all scenarios by convolving the slit profile scaled to each value with the flat field
using the pre-existing slit_flat_convolve function of the ghost class.

It then loops through every order and section taking the maximum of the convolution and corresponding
scale for each section, in order to determine what is the best scale for each case.

The resulting measured for all sections are then fitted in a least squares approach
to determine the best model for the spectrograph. This fit is weighted by the maximum convolution value
to ensure the edges of the chip where the flat flux is low are down weighted.
