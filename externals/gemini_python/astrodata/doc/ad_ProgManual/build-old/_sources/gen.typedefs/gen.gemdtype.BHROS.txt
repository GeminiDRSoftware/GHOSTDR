
BHROS Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    BHROS

Source Location 
    ADCONFIG_Gemini/classifications/types/BHROS/gemdtype.BHROS.py

.. code-block:: python
    :linenos:

    class BHROS(DataClassification):
        name="BHROS"
        usage = ""
        typeReqs= []
        phuReqs= {}
        parent = "GEMINI"
        requirement = PHU(INSTRUME='bHROS')
    
    newtypes.append(BHROS())



