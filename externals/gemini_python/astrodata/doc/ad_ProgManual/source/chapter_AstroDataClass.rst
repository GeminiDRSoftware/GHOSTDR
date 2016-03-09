AstroData Class Reference
!!!!!!!!!!!!!!!!!!!!!!!!!

.. toctree::
    
The following is information about the ``AstroData`` class. For descriptions of
arguments shown for the class constructor, see ``AstroData.__init__(..)``.  This
documentation is generated in part from in-source docstrings. 

To import the ``AstroData`` class use::

    from astrodata import AstroData


AstroData Class
@@@@@@@@@@@@@@@

.. autoclass:: astrodata.AstroData

Basic Functions
@@@@@@@@@@@@@@@

AstroData Constructor
#####################

.. toctree::
    
.. automethod:: astrodata.AstroData.__init__

append(..)
##########

.. toctree::

.. automethod:: astrodata.AstroData.append

close(..)
#########

.. toctree::

.. automethod:: astrodata.AstroData.close

insert(..)
##########

.. toctree::

.. automethod:: astrodata.AstroData.insert

info(..)
########

.. toctree::

.. automethod:: astrodata.AstroData.info

write(..)
#########

.. toctree::

.. automethod:: astrodata.AstroData.write

Type Information
@@@@@@@@@@@@@@@@

.. toctree::

.. automethod:: astrodata.AstroData.type

.. automethod:: astrodata.AstroData.status


Header Manipulations
@@@@@@@@@@@@@@@@@@@@

.. toctree::

Manipulations of headers, specifically retrieving and setting key-value
pair settings in the header section of header-data units
can be done directly using the AstroData header manipulation functions
which cover both PHU and extension headers.
For higher level metadata which is available for all types in the tree
in a properly constructed configuration space, the metadata is retrieved with 
descriptor functions, accessed as members of the AstroData object.

To retrieve or set meta-data not covered by descriptors, one must
read and write key-value pairs to the HDU headers at the lower-level. AstroData
offers three pairs of functions for getting and setting header values, for each
of three distinct cases.  While it is possible to use the pyfits.Header directly
(available via "ad[..].header"), it is preferrable to use the AstroData calls
which allow AstroData to keep type information up to date, as well as to update
any other characteristics of the AstroData object which may need to be
maintained when the dataset is changed.

The three distinct pairs of header access functions serve the following
purposes:

+ set/get headers in PHU.
+ set/get headers in the single extension of a "single-HDU AstroData
  object".
+ set/get headers in an extension of a multi-HDU (aka "multi-extension") 
  AstroData instance. This requires specifying the extension index, and
  cannot be used to modify the PHU. HDU #0 is the first real 
  header-data section in the MEF.

Set/Get PHU Headers
###################

.. toctree::

.. automethod:: astrodata.AstroData.phu_get_key_value

.. automethod:: astrodata.AstroData.phu_set_key_value

Set/Get Single-HDU Headers
##########################

.. toctree::

.. automethod:: astrodata.AstroData.get_key_value

.. automethod:: astrodata.AstroData.set_key_value


Iteration and Subdata
@@@@@@@@@@@@@@@@@@@@@

.. toctree::
   
Overview
########
 
.. toctree::
    
    gen.ADMANUAL-ADSubdata.rst
    
count_exts(..)
##############

.. toctree::

.. automethod:: astrodata.AstroData.count_exts


The [] Operator
###############

.. toctree::

.. automethod:: astrodata.AstroData.__getitem__

data attribute
##############
 
.. toctree::
     
.. autoattribute:: astrodata.AstroData.data

descriptors attribute
#####################
 
.. toctree::
     
.. autoattribute:: astrodata.AstroData.descriptors

filename attribute
##################
 
.. toctree::
     
.. autoattribute:: astrodata.AstroData.filename

header attribute
################

.. toctree::

.. autoattribute:: astrodata.AstroData.header

headers attribute
#################

.. toctree::

.. autoattribute:: astrodata.AstroData.headers

hdulist attribute
#################
 
.. toctree::
     
.. autoattribute:: astrodata.AstroData.hdulist

phu attribute
#############
 
.. toctree::
     
.. autoattribute:: astrodata.AstroData.phu


Renaming an Extension
#####################

.. toctree::

.. automethod:: astrodata.AstroData.rename_ext
