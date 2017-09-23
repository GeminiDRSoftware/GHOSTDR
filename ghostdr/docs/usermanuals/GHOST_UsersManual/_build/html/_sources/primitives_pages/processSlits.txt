.. primitive1:

.. processSlits:

processSlits
============================

Purpose
-------

This primitive replaces CR-affected pixels in each individual slit viewer image
(taken from the current stream) with their equivalents from the median frame of
those images, and extracts slit object profiles to flux-weight the images'
offsets (relative to the coincident science frame) from which it computes the
mean exposure epoch.

Inputs and Outputs
------------------

``flat``: string or None (default: None) - The name of the processed slit viewer
flat field image.  If not provided, or if set to None, the flat will be taken
from the ``processed_slitflat`` override_cal command line argument.

``slit``: string or None (default: None) - The name of the bias-/dark-corrected
median slit viewer frame.  If not provided, or if set to None, the frame will be
taken from the ``slitStream`` argument.

``slitStream``: string or None (default: None) - The name of the stream from
which to access the bias-/dark-corrected median slit viewer frame (only the
first AstroData instance will be used).  If not provided, or if set to None, the
frame will be taken from the ``processed_slit`` override_cal command line
argument.

Algorithm
---------

The incoming stream of AstroData objects, representing
bias-/dark-corrected slit viewer frames, is first used to produce a per-pixel
median absolute deviation (MAD) frame for clipping cosmic ray (CR) outliers.
The median frame is subtracted from each individual frame in turn, and where the
residual exceeds 20x the MAD, a CR is detected.  Each CR-affected pixel is
replaced with its equivalent from the median slit viewer frame.

Then, for each CR-corrected individual frame, a ``SlitViewer`` object
is instantiated for each spectrograph arm using the frame and the passed-in flat
field image.  The ``SlitViewer`` objects then each provide, via the
``object_slit_profiles()`` method, a 2 element array where each element is a
1-dimensional array representing a single object profile, un-normalised,
sky-corrected, and summed across (not along) the slit.

In the case of high resolution data, the Th/Ar profiles are discarded, after
which the flux for the remaining object profiles, for both arms, is summed
together.  This individual frame sum is then scaled by its percentage overlap
with the science frame (computed from UTC start/stop times for both the
individual frame and the science frame), which in turn is used to weight the
frame's offset (in seconds) from the start of the science exposure.  [Note that
only the science frame's UTC start/stop times are present for this process, not
its frame data.]

Finally the sum of the weighted offsets is divided by the sum of the weights to
provide the mean offset from the start of the science exposure.  This is added
to the science exposure UTC start time to produce the mean exposure epoch, which
in turn is written into the header of each CR-corrected output frame, in the
``AVGEPOCH`` keyword.

Issues and Limitations
----------------------

For short duration science exposures there may be too few coincident slit frames
to produce a meaningful MAD frame for CR detection.
