.. primitive1:

.. promote:

promote
============================

Purpose
-------
This primitive is responsible for taking a slit viewer MEF (that results as the
output from the observation bundle splitter primitive) and "promoting" all the
slit viewer extensions therein to full AstroData objects with various
manipulations so that they will pass various checks and otherwise work nicely
with the rest of the recipe.

Note that this primitive is only meant to be run on raw (unprocessed) slit
viewer MEFs produced by the splitter primitive (or by the simulator with the
-split option supplied).

Inputs and Outputs
------------------

``promote`` takes no particular configuration inputs.

Algorithm
---------

The ``promote`` primitive performs these manipulations on the slit viewer
extensions:

- It uses AstroData slicing (also known as "sub-data") to initially produce a
  separate AstroData object for each extension, then also makes a copy of the
  PHU for each resultant AstroData object (since sub-data references a single
  copy of the original PHU by default).  In this manner, fully independent
  AstroData objects are created representing each extension.  This is done to
  bypass the extension count check we inherit from
  ``GMOSPrimitives.validateData``, which is suitable for processing GHOST's red
  and blue arm data, but not for the slit viewer MEFs which may contain any
  number of extensions.
- The above manipulation causes the resulting AstroData objects to all have the
  same ORIGNAME in their PHUs.  This becomes an issue when we leverage the
  existing ``StackPrimitives.stackFrames`` primitive, as it must first (being
  IRAF-based) write its inputs to disk.  As it uses ORIGNAME to do so, clashes
  would result, so we make ORIGNAME unique by simply appending an incrementing
  number.
- Finally, we also reset the EXTVER of all extensions in the AstroData objects
  resulting from the first manipulation back to 1 (since now there is only a
  single extension in each AstroData object).

Issues and Limitations
----------------------

None.
