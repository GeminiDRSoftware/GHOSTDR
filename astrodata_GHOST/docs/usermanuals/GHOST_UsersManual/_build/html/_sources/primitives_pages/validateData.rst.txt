.. primitive1:

.. validateData:

validateData
============================

Purpose
-------

This primitive is responsible for ensuring the data passed is GHOST data and has
correctly formatted keywords.  (We use the "prepare" superprimitive provided by
Gemini which invokes validateData before standardizeHeaders and
standardizeStructure, so validateData must not assert things set by either of
those.)

Inputs and Outputs
------------------

``validateData`` takes no particular configuration inputs.

Algorithm
---------

None.

Issues and Limitations
----------------------

At the moment, our version of this primitive does nothing other than write a
timestamp into the headers.  In future we may want to check that the right
number of extensions is present (but only for normal science detector frames as
slit viewer frames can have any number), and/or confirming that only certain
binning shape is used.  These are the sorts of assertions other instruments
have used validateData for.
