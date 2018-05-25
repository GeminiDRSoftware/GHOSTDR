:tocdepth: 2

.. primitives:

*****************
Primitive Classes
*****************

.. autoclass:: ghostdr.ghost.primitives_ghost.GHOST
   :show-inheritance:

.. autoclass:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect
   :show-inheritance:

.. autoclass:: ghostdr.ghost.primitives_ghost_bundle.GHOSTBundle
   :show-inheritance:

.. autoclass:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit
   :show-inheritance:

***************
Main Primitives
***************

``addWavelengthSolution``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.addWavelengthSolution

Inputs and Outputs
------------------

Each extension of the input AstroData object has a ``.WAVL`` data extension
added to it. This extension gives wavelength values in a one-to-one
correspondence with the data dispersion axis.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

- Wavelength units are currently not put into the output header. This is to
  be rectified.

``applyFlatBPM``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.applyFlatBPM

Inputs and Outputs
------------------

This primitive does not change the input AstroData structure..

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If the science frame and flat field BPMs have a different shape and/or
extension structure, the primitive will kill the RecipeSystem rather than
attempt to figure out what has gone wrong.

``barycentricCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.barycentricCorrect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``clipSigmaBPM``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.clipSigmaBPM

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``darkCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.darkCorrect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``extractProfile``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.extractProfile

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``interpolateAndCombine``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.interpolateAndCombine

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``findApertures``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.findApertures

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``flatCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.flatCorrect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``formatOutput``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.formatOutput

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``responseCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.responseCorrect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``splitBundle``
============================

:class:`ghostdr.ghost.primitives_ghost_bundle.GHOSTBundle`

.. automethod:: ghostdr.ghost.primitives_ghost_bundle.GHOSTBundle.splitBundle

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``standardizeStructure``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.standardizeStructure

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``tileArrays``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.tileArrays

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

**********************
Calibration Primitives
**********************

``CRCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`

.. automethod:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit.CRCorrect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``fitWavelength``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.fitWavelength

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``processSlits``
============================

:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`

.. automethod:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit.processSlits

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``stackFrames``
============================

:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`

.. automethod:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit.stackFrames

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

****************
Helper Functions
****************

These functions reside in the main primitive modules, but are not intended to
be called as primitives by end users.

``_compute_barycentric_correction``
===================================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect._compute_barycentric_correction

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_get_polyfit_filename``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect._get_polyfit_filename

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_interp_spect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect._interp_spect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_mad``
============================

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autofunction:: _mad

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_rebin_ghost_ad``
============================

:class:`ghostdr.ghost.primitives_ghost.GHOST`

.. automethod:: ghostdr.ghost.primitives_ghost.GHOST._rebin_ghost_ad

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_regrid_spect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect._regrid_spect

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_request_bracket_arc``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect._request_bracket_arc

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_total_obj_flux``
============================

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autofunction:: _total_obj_flux

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.

``_write_newfile``
============================

.. currentmodule:: ghostdr.ghost.primitives_ghost_bundle
.. autofunction:: _write_newfile

.. automodule:: ghostdr.ghost.primitives_ghost_bundle
   :members: _write_newfile

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

AstroData Tag(s)
-----------------

List the types/tags (or sets of tags)

Algorithms
----------

Algorithm details here (if applicable)

Issues and Limitations
----------------------

If applicable.
