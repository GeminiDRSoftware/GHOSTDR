
NIRI_SPECT Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    NIRI_SPECT

Source Location 
    ADCONFIG_Gemini/classifications/types/NIRI/gemdtype.NIRI_SPECT.py

.. code-block:: python
    :linenos:

    
    class NIRI_SPECT(DataClassification):
        name="NIRI_SPECT"
        usage = "Applies to any spectra from the NIRI instrument."
        parent = "NIRI"
        requirement = AND([ISCLASS('NIRI'),
                           PHU(FILTER3='(.*?)grism(.*?)'),
                           NOT(ISCLASS("NIRI_DARK"))
                           ])
                           
    
    newtypes.append(NIRI_SPECT())



