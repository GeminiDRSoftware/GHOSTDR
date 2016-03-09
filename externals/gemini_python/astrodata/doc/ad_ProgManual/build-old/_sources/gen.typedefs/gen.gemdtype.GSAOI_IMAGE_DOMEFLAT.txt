
GSAOI_IMAGE_DOMEFLAT Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    GSAOI_IMAGE_DOMEFLAT

Source Location 
    ADCONFIG_Gemini/classifications/types/GSAOI/gemdtype.GSAOI_IMAGE_DOMEFLAT.py

.. code-block:: python
    :linenos:

    class GSAOI_IMAGE_DOMEFLAT(DataClassification):
        name="GSAOI_IMAGE_DOMEFLAT"
        usage = """
            Applies to all imaging domeflat flat datasets from the GSAOI
            instrument
            """
        parent = "GSAOI_IMAGE_FLAT"
        requirement = AND([  ISCLASS("GSAOI_IMAGE_FLAT"),
                             PHU(OBJECT="(([Dd]omeflat)+?)")  ])
    
    newtypes.append(GSAOI_IMAGE_DOMEFLAT())



