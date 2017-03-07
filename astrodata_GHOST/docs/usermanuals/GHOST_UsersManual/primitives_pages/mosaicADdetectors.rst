.. primitive1:

.. mosaicADdetectors:

mosaicADdetectors
============================

Purpose
-------

This primitive mosaics the SCI frames of the input images, along with the VAR
and DQ frames if they exist.

Inputs and Outputs
------------------

``tile``: bool (default: False), tile images instead of mosaic (mosaic includes
transformations)

``dq_planes``: bool (default: False), transform the DQ image, bit plane by bit
plane

Algorithm
---------

This primitive overrides the one of the same name inherited from
``GEMINIPrimitives``.  The override logic is identical to the original except
for two hacks introduced just prior to calling the mosaicing routine,
``gemini_mosaic_function``, and undone immediately after.  These are necessary
in order to trick the underlying mosaicing machinery into working for GHOST
data, and are as follows:

- First, we add 'GSAOI' to the internal array of types for our input AstroData
  object (which contains the extensions we want to mosaic).  This bypasses the
  GSAOI-/GMOS-only check built into ``gemini_mosaic_function``.  Afterwards we
  reset the type list back to normal by invoking ``AstroData.refresh_types`` on
  the resulting, mosaiced, AstroData object.
- Secondly, we set the ``INSTRUME`` PHU keyword to
  ``../../../astrodata_GHOST/ADCONFIG_GHOST/lookups``.  This tricks
  ``set_geo_values`` (called by ``gemini_mosaic_function``) into referencing the
  ``geometry_conf.py`` file under our own ``astrodata_GHOST/lookups`` folder. As
  soon as the mosaicing routine returns we reset ``INSTRUME`` back to 'GHOST'
  (and in fact this must be done before the call to ``refresh_types`` mentioned
  above).

Until such time as ``astrodata_GHOST`` has been merged into ``astrodata_Gemini``
(like GMOS and GSAOI have been), the above two hacks will be necessary in order
for GHOST to make use of Gemini's mosaicing mechanism.

Issues and Limitations
----------------------

Requires ugly hacks until integration with ``astrodata_Gemini`` is completed.
