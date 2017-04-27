.. primitive1:

.. switchTo:

switchTo
============================

Purpose
-------
This primitive makes the named stream the current stream, such that all
subsequently invoked primitives will take as input the AstroData objects in that
stream, as well as write their output AstroData objects to that stream.

Inputs and Outputs
------------------

The only parameter ``switchTo`` takes is ``streamName``, the name of the stream
to which to "shift" execution.  Equivalently, the name of the stream which to
make current.

Issues and Limitations
----------------------

Be careful not to specify the parameter ``streamName`` as simply ``stream`` as
this is a reserved parameter name with special meaning to the recipe framework
and odd behaviour results when used.
