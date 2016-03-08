#test
from astrodata.mkro import *


def name_of_primitive(*args, **argv):
    """
    First argument is an astrodata object.  Then any other arguments
    that is obtained by the primitive through the ReductionContext, "rc",
    can be passed.

    Example:
       from pifghost import name_of_primitive
       name_of_primitive(ad, threshold=50.)
    """

    ro = mkRO(astrotype="GHOST_SPECT", copy_input=True, 
              args=args, argv=argv)
    ro.runstep("nameOfPrimitive", ro.context)
    outputs = ro.context.get_outputs(style="AD")
    if len(outputs)==0:
        return None
    elif len(outputs)==1:
        return outputs[0]
    else:
        return outputs
    
    
