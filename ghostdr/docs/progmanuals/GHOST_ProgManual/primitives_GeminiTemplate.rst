:tocdepth: 2

.. primitives:

***************
Main Primitives
***************

applyFlatBPM
============================

.. autoclass:: ghostdr.ghost.primitives_ghost_spect.GHOSTSpect
   :members: applyFlatBPM
   :show-inheritance:

Purpose
-------

The profile extraction routine requires that the science frame apertures be
extracted using the aperture model determined from the flat field. However,
as the flat field corrected has not been applied to the science frame by this
point in the reduction sequence, pixels which are marked bad in the flat field
may be included in the science frame aperture extraction.

This primitive is designed to negate that problem by folding the flat field
bad pixel mask (BPM) into the science frame BPM before any aperture extraction
or flat field correction is performed.

Inputs and Outputs
------------------

Describe state of file before and after processing, including which types/tags
will be found/applied before and after

Input parameters
----------------

Table of parameter name - meaning

AstroData Type(s)
-----------------

Just list the types/tags (or sets of tags)

Inheritance and Primitive Set
-----------------------------

List primitive classes that this class inherits from

Location
--------

Import path

.. class:: astrodata_GHOST.RECIPES_GHOST.primitives.primitives_GHOST.GHOSTPrimitives

Algorithms
----------

This primitive uses a standard ``numpy.bitwise_or`` call to combine the BPMs.

Issues and Limitations
----------------------

If the science frame and flat field BPMs have a different shape and/or
extension structure, the primitive will kill the RecipeSystem rather than
attempt to figure out what has gone wrong.

**********************
Calibration Primitives
**********************

Primitive #2
============

