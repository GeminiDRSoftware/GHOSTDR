#
#                                                                  gemini_python
#
#                                                                astrodata.utils
#                                                                  debugmodes.py
# ------------------------------------------------------------------------------
# $Id: debugmodes.py 5274 2015-06-11 14:39:37Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5274 $'[11:-2]
__version_date__ = '$Date: 2015-06-12 00:39:37 +1000 (Fri, 12 Jun 2015) $'[7:-2]
# ------------------------------------------------------------------------------
_throw_descriptor_exception = False

def set_descriptor_throw(flag = True):
    """ This function stores a flag the descriptor system uses to throw exceptions
    when descriptors do, as opposed to the default behavior of catching the exception
    and returning None. """
    global _throw_descriptor_exception
    _throw_descriptor_exception = flag
     
def get_descriptor_throw():
    global _throw_descriptor_exception
    return _throw_descriptor_exception
