.. primitive1:

.. correctSlitCosmics:

correctSlitCosmics
============================

Purpose
-------

This primitive replaces CR-affected pixels in each individual slit viewer image
(taken from the current stream) with their equivalents from the median frame of
those images.

Inputs and Outputs
------------------

``correctSlitCosmics`` takes no particular configuration inputs.

Algorithm
---------

The incoming stream of AstroData objects, representing
bias-/dark-corrected slit viewer frames, is first used to produce a per-pixel
median absolute deviation (MAD) frame for clipping cosmic ray (CR) outliers.
A median frame is then separately produced (from the inputs) and subtracted
from each in turn.  CR-impacted pixels are those where the difference residual
exceeds 20x the MAD, which are then replaced with their equivalents from the
median slit viewer frame.

Issues and Limitations
----------------------

For short duration science exposures there may be too few coincident slit frames
to produce a meaningful MAD frame for CR detection.
