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

This primitive does not change the input AstroData structure.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

None known.

``barycentricCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.barycentricCorrect

Inputs and Outputs
------------------

This primitive does not change the input AstroData structure.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

None known.

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

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

None known.

``darkCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.darkCorrect

Inputs and Outputs
------------------

No changes are made to the structure of the input AstroData object.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

As mentioned above, there must be a one-to-one correspondence between the input
AstroData objects and the list of dark files to be applied. There is currently
no way to broadcast a single dark to multiple AstroData inputs.

``extractProfile``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.extractProfile

Inputs and Outputs
------------------

After processing with this primitive, the AstroData object changes from being
an image file to a file containing extracted flux values (i.e., uncalibrated
spectra). The original image data is *not* retained, although it can be kept
by forcing a ``write_result`` on the last primitive invoked before calling
``extractProfile``.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Algorithms
----------

The concept behind :any:`polyfit <polyfit>` is described in :ref:`polyfit-core`.

Issues and Limitations
----------------------

None known.

``interpolateAndCombine``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.interpolateAndCombine

Inputs and Outputs
------------------

The un-interpolated spectrum data is left in place. The interpolated
data is provided as an additional extension to the end of the input
AstroData file.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

The only implemented wavelength scale is ``'loglinear'``. A ``'linear'`` option
should be provided at a later date.

``findApertures``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.findApertures

Inputs and Outputs
------------------

Each extension of the input AstroData gains a new attribute, ``.XMOD``, which
describes the polyfit model which represent the aperture locations.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Algorithms
----------

The algorithm behind generating a polyfit model is described in
:ref:`polyfit-models`.

Issues and Limitations
----------------------

This primitive currently contains what is described as an 'attempt' to remove
the worst cosmic rays from teh data before it is passed along to the model
fitting routines. This is done by applying a :math:`5\times 5` median filter
to the data, and then replacing any points which have percetage difference
between the real data and median filter value with the median filter value.
Furthermore, only data points above the average data value are replaced.
Presumably, this will be removed when an improved extraction routine is
developed, which should remove cosmic rays intrinsically.

``flatCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.flatCorrect

Inputs and Outputs
------------------

This primitive does not change the structure of the AstroData input(s).

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

.. warning::
            While the primitive is working, it has been found that the
            underlying algorithm is flawed. A new algorithm is being developed.

``formatOutput``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.formatOutput

Inputs and Outputs
------------------

Additional extensions will be added to the input AstroData objects, one
extension per requested data product.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

The options for additional data products are deliberately sequential, i.e.,
if the user asks for the last additional data product in the reduction
sequence, they will also get all the previous additional data products. This is
to guarantee that the same data product will always occur at the same
extension number in output files (given that extensions can't be explicitly
labelled).

``responseCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.responseCorrect

Inputs and Outputs
------------------

This primitive does not alter the structure of the input AstroData object(s).

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

This primitive will need to be updated once the Gemini calibration service
for flux standards becomes available.

``splitBundle``
============================

:class:`ghostdr.ghost.primitives_ghost_bundle.GHOSTBundle`

.. automethod:: ghostdr.ghost.primitives_ghost_bundle.GHOSTBundle.splitBundle

Inputs and Outputs
------------------

GHOST observations come packaged in a 'bundle' FITS file, containing the
exposure data from all instrument cameras during a given exposure. This
primitive breaks that bundle file into its core components:

- The red camera data;
- The blue camera data;
- Any slit viewer exposure(s) taken during the main observation.

Each output uses the original bundle file name, modified to identify which
camera it corresponds to (as well as any filename suffix passed).

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_bundle
.. autoattribute:: GHOSTBundle.tagset

Algorithms
----------

This primitive makes use of the :func:`_write_newfile <_write_newfile>` helper
function.

Issues and Limitations
----------------------

None known.

``standardizeStructure``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.standardizeStructure

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

None known.

``tileArrays``
============================


:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.tileArrays

Inputs and Outputs
------------------

The four individual image extensions of the input AstroData object(s) are
stitched together into a single image frame.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_spect
.. autoattribute:: GHOSTSpect.tagset

Issues and Limitations
----------------------

This primitive contains an internal hepler function, ``simple_mosaic_function``,
to do the actual tiling. This function will (probably) be placed somewhere in
the higher-level Gemini primitives at a later date.

**********************
Calibration Primitives
**********************

``CRCorrect``
============================

:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`

.. automethod:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit.CRCorrect

Inputs and Outputs
------------------

The structure of the input slit images is not changed.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autoattribute:: GHOSTSlit.tagset

Issues and Limitations
----------------------

None known.

``fitWavelength``
============================

:class:`ghostdr.ghost.primitives_ghost_spect.GHOSTSpect`

.. automethod:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect.fitWavelength

Inputs and Outputs
------------------

The arc file data extension(s) will gain a new attribute, ``.WFIT``, which
stores the wavelength solution.

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autoattribute:: GHOSTSlit.tagset

Algorithms
----------

.. note::
   This needs to point to the ``read_lines_and_fit`` documentation when it
   exists.

Issues and Limitations
----------------------

None known.

``processSlits``
============================

:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`

.. automethod:: ghostdr.ghost.primitives_ghost_slit.GHOSTSlit.processSlits

Inputs and Outputs
------------------

This primitve does not change the structure of the input AstroData object(s).

AstroData Tag(s)
-----------------

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autoattribute:: GHOSTSlit.tagset

Issues and Limitations
----------------------

None known.

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

.. currentmodule:: ghostdr.ghost.primitives_ghost_slit
.. autoattribute:: GHOSTSlit.tagset

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
