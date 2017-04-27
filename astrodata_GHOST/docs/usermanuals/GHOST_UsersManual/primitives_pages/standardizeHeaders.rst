.. primitive1:

.. standardizeHeaders:

standardizeHeaders
============================

Purpose
-------
This primitive implements part of the necessary interface required by the
``StandardizePrimitives.prepare`` "super-primitive" inherited from Gemini which,
among other things, we use for its calls to ``validateData`` and
``markAsPrepared``.  Gemini defines both of those so we don't have to, but not
``standardizeHeaders``; however, it does define ``standardizeGeminiHeaders``
which works fine for us so, ironically, we just make this primitive call that!

Inputs and Outputs
------------------

The usual stream of AstroData objects which, upon output, have had their headers
standardized for trouble-free processing through the GHOST pipeline recipes.

Issues and Limitations
----------------------

None.
