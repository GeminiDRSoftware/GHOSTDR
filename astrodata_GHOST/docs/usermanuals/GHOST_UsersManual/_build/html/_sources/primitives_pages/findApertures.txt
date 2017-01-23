.. primitive1:

.. rejectCosmicRays:

findApertures
============================

Purpose
-------
This primitive will take an existing ``polyfit`` model for the aperture
locations in a frame, and fit the model to the observed data.

Note that this primitive is only run on prepared flat fields (i.e. images
that have been run through the ``makeProcessedFlatG`` recipe). The primitive
will abort if passed any other type of file.

Inputs and Outputs
------------------

``findApertures`` takes no particular configuration inputs.

Algorithm
---------

The RecipeSystem will automatically determine the appropriate initial polyfit
model, based on the arm, resolution and date of the observation.

.. note:: Date selection of polyfit models has yet to be implemented.

The primitive will instantiate an ``Arm`` object from the ``polyfit.ghost``
module, which is included as a part of ``astrodata_GHOST``. This ``Arm``
contains all the functions for fitting the aperture positions of the data.

The ``Arm`` object then makes the following three steps:

- An initial model of the spectrograph is constructed based on the parameters
  read-on from the lookup system;
- The reduced flat-field is convolved with the slit profile of the instrument,
  determined from the slit viewing camera
- The initial model is fitted to the result of the convolved flat field, which
  represents the location of the middle of each order. The result of this fit
  is then used in the spectra extraction as the basis for the location of the
  orders.

A ``polyfit`` model FITS file is then written out to the calibration system.


Issues and Limitations
----------------------

There is a placeholder retrieval function for getting back a fitted
``polyfit`` model, but there are no primitives yet for applying this model
to data.

Future versions of this fitting procedure will include an extra polynomial
describing the change in the spatial scale of the fiber images as a function
of order location on the CCD. The current lack of such feature currently
results in a slight innaccuracy of the fit on the edge of orders which are
close together. This is due to the fact that the convolution with a fixed profile
results is slight overlap at those locations.

