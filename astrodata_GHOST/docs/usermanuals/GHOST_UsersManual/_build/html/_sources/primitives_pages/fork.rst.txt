.. primitive1:

.. fork:

fork
============================

Purpose
-------
This primitive creates a new stream by copying the current stream's inputs to
the outputs of the new stream.  Has the same effect as (but without the disk
write penalty incurred by) the following construct::

	addToList(purpose=save_to_disk)
	getList(purpose=save_to_disk, to_stream=new_stream_name)

Inputs and Outputs
------------------

The only parameter ``fork`` takes is ``newStream``, the name of the new stream
to be formed.

Issues and Limitations
----------------------

Be careful not to specify the parameter ``newStream`` as simply ``stream`` as
this is a reserved parameter name with special meaning to the recipe framework
and odd behaviour results when used.
