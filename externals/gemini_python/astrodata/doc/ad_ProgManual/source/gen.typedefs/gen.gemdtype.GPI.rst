
GPI Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    GPI

Source Location 
    ADCONFIG_Gemini/classifications/types/GPI/gemdtype.GPI.py

.. code-block:: python
    :linenos:

    
    class GPI(DataClassification):
        name="GPI"
        usage = "Applies to datasets from the GPI instrument"
        parent = "GEMINI"
        requirement = PHU(INSTRUME='GPI')
    
    newtypes.append(GPI())



