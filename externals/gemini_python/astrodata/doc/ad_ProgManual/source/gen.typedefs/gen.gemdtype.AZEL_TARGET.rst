
AZEL_TARGET Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    AZEL_TARGET

Source Location 
    ADCONFIG_Gemini/classifications/types/gemdtype.AZEL_TARGET.py

.. code-block:: python
    :linenos:

    
    class AZEL_TARGET(DataClassification):
        name="AZEL_TARGET"
        usage = "Data taken on a target in the AZEL_TOPO co-ordinate system. This means that the telescope pointing position was defined as an (Azimuth, Elevation) position, not an (RA, Dec) position. The telesope is not tracking a celestial source position. This is mostly used for daytime engineering tests."
        
        parent = "GEMINI"
        requirement = PHU(FRAME='AZEL_TOPO') 
    
    newtypes.append(AZEL_TARGET())



