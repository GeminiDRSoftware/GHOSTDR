.. primitive1:

.. processSlits:

processSlits
============================

Purpose
-------

For each CR-corrected slit viewer frame taken from the current stream, this
primitive extracts the object profiles within and uses them to flux-weight the
image's offset (relative to the coincident science frame) and then uses that to
compute the mean exposure epoch for the entire series of inputs.

Inputs and Outputs
------------------

``flat``: string or None (default: None) - The name of the processed slit viewer
flat field image.  If not provided, or if set to None, the flat will be taken
from the ``processed_slitflat`` override_cal command line argument.

Algorithm
---------

For each individual CR-corrected frame taken from the input, a ``SlitViewer``
object is instantiated for each spectrograph arm using the frame and the passed-
in flat field image.  The ``SlitViewer`` objects then each provide, via the
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

The UTC start/stop times of the coincident science frame are currently obtained
from the UTSTART/UTEND PHU keywords, having been written there by our splitter
primitive (or the simulator with the -split option used).  In the context of a
standalone slit viewer frame, this may be considered improper usage of these
keywords.
