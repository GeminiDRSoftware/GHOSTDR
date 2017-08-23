:tocdepth: 2

.. primitives:

**********
Primitives
**********

Core Primitives
===============

Sub-header
----------

Sub-subheader
+++++++++++++

Calibration Primitives
======================

These primitives are inherited by GHOSTPrimitives.

.. note:: Most of these primitives are currently undocumented. This will
          be rectified once the Gemini code refactor is completed.

.. module:: astrodata_GHOST.RECIPES_GHOST.primitives.primitives_GHOST_calibration
.. autoclass:: GHOST_CalibrationPrimitives
    :members:
    :private-members:

**********
Primitives
**********

applyFlatBPM
============================

.. autoclass:: astrodata_GHOST.RECIPES_GHOST.primitives.primitives_GHOST.GHOSTPrimitives
   :members: applyFlatBPM

Purpose
-------

The profile extraction routine requires that the profile be extracted using
the trace determine from the flat field, but *before* flat field correction
is applied. This means that the science frame will have no knowledge of bad
pixels in the flat used to establish the trace.

This primitive is designed to negate that problem by folding the flat field
bad pixel mask (BPM) into the science frame BPM before any aperture extraction
or flat field correction is performed.

Inputs and Outputs
------------------

Input parameters
----------------

AstroData Type(s)
-----------------

Inheritance and Primitive Set
-----------------------------

Location
--------

Algorithms
----------

Issues and Limitations
----------------------

Primitive #2
============

