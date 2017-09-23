.. primitive1:

.. standardizeStructure:

standardizeStructure
============================

Purpose
-------
This primitive is responsible for massaging input data frames into a format
compatible with the remainder of the downstream primitives.

When passed a series of slit viewer MEFs (such as results from the observation
bundle splitter primitive, or from the simulator with the -split option sup-
plied), it "promotes" all the slit viewer extensions therein to full AstroData
objects with various manipulations so that they will pass various checks and
otherwise work nicely with the rest of the recipe.

When passed other input (such as normal science detector frames) nothing is done
except the addition of the obligatory timestamp keyword.

Inputs and Outputs
------------------

``standardizeStructure`` takes no particular configuration inputs.

Algorithm
---------

For slit viewer data, the ``standardizeStructure`` primitive performs these man-
ipulations:

- It uses AstroData slicing (also known as "sub-data") to initially produce a
  separate AstroData object for each extension, then also makes a copy of the
  PHU for each resultant AstroData object (since sub-data references a single
  copy of the original PHU by default).  In this manner, fully independent
  AstroData objects are created representing each extension.
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
