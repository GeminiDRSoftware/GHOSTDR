ReductionContext Class Reference
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The following is information about the ReductionContext class. When writing
primitives the reduction context is passed into the primitive as the sole
argument (generally named ``rc`` by
Gemini conventions and in addition to the ``self`` argument).
This object is used by the primitive to both get inputs
and store outputs, as well as to communicate with subsystems
like the calibration queries system or list keeping for stacking.

Parameter and Dictionary Features
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

The "in" operator: contains(..)
###############################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.__contains__

Dataset Streams: Input and Output Datasets
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

get_inputs(..)
##############

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_inputs

get_inputs_as_astrodata(..)
###########################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_inputs_as_astrodata

get_inputs_as_filenames(..)
###########################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_inputs_as_filenames

get_stream(..)
##############

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_stream

get_reference_image(..)
#######################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_reference_image

report_output(..)
#################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.report_output

switch_stream(..)
#################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.switch_stream


Calibrations
@@@@@@@@@@@@

get_cal(..)
###########

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_cal

rq_cal(..)
##########

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.rq_cal

Stacking
@@@@@@@@

rq_stack_get(..)
################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.rq_stack_get

rq_stack_update(..)
###################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.rq_stack_update

Lists
@@@@@

list_append(..)
###############

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.list_append

get_list(..)
############

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.get_list


Utility
@@@@@@@

prepend_names(..)
#################

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.prepend_names


run(..)
#######

.. toctree::

.. automethod:: recipe_system.reduction.reductionContext.ReductionContext.run
