
GCAL_IR_ON Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    GCAL_IR_ON

Source Location 
    ADCONFIG_Gemini/classifications/types/gemdtype.GCAL_IR_ON.py

.. code-block:: python
    :linenos:

    class GCAL_IR_ON(DataClassification):
        name="GCAL_IR_ON"
        usage = "Indicates that the GCAL IR flat field lamp is on and the shutter is open. This is typically referred to as a lamp-on flat"
        
        parent = "GEMINI"
        requirement = PHU(GCALLAMP='IRhigh') & PHU(GCALSHUT='OPEN')
    
    newtypes.append(GCAL_IR_ON())



