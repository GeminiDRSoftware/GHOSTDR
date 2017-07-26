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
