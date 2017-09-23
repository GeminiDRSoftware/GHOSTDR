.. primitive1:

.. stackFrames:

stackFrames
============================

Purpose
-------
This primitive is simply a wrapper around the standard
``StackPrimitives.stackFrames`` primitive which exposes extra parameters to the
underlying ``gemcombine.cl`` IRAF script, and allows increasing the IRAF
verbosity.

Inputs and Outputs
------------------

``hsigma``: float or None (default: None) - expose the hsigma parameter of the
underlying IRAF gemcombine script, allowing it to be set within a recipe

``hthresh``: float or None (default: None) - expose the hthreshold parameter of
the underlying IRAF gemcombine script, allowing it to be set within a recipe

``lsigma``: float or None (default: None) - expose the lsigma parameter of the
underlying IRAF gemcombine script, allowing it to be set within a recipe

``lthresh``: float or None (default: None) - expose the lthreshold parameter of
the underlying IRAF gemcombine script, allowing it to be set within a recipe

``pclip``: float or None (default: None) - expose the pclip parameter of the
underlying IRAF gemcombine script, allowing it to be set within a recipe

``sigscale``: float or None (default: None) - expose the sigscale parameter of
the underlying IRAF gemcombine script, allowing it to be set within a recipe

``verbose``: <any> or None (default: None) - set the level of iraf verbosity

Issues and Limitations
----------------------

The default value for ``reject_method`` (a parameter exposed by the underlying
``StackPrimitives.stackFrames`` primitive) is 'avsigclip' (inherited from
Gemini) which is probably not what we want in most cases.  If necessary (or more
convenient), we could override the default in this wrapper, or in the
``parameters_GHOST.py`` file.
