#
#                                                                  gemini_python
#
#                                                                astrodata.utils
#                                                                   adinspect.py
# ------------------------------------------------------------------------------
# $Id: adinspect.py 5274 2015-06-11 14:39:37Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5274 $'[11:-2]
__version_date__ = '$Date: 2015-06-11 04:39:37 -1000 (Thu, 11 Jun 2015) $'[7:-2]
# ------------------------------------------------------------------------------
import inspect

def get_primitives(primsetInstance):
    def pred(tob):
        return inspect.isgeneratorfunction(tob)
    membs = inspect.getmembers(primsetInstance, pred)
    names = [ a[0] for a in membs]
    return names

def get_descriptors(classname=None, astrotype=None):
    from ..interface.Descriptors import centralCalculatorIndex
        
    classname = centralCalculatorIndex[astrotype]
    cnp = classname.split(".")
    exec("import "+cnp[0])    
    calcobj = eval(classname)
    membs = inspect.getmembers(calcobj,inspect.ismethod)
    ret = {}
    for memb in membs:
        newdict = {}
        print "26",memb[0],repr(type(memb[1]))
        if memb[0][0] == "_":
            continue
        ret.update({memb[0]:newdict})
        newdict.update({"name":memb[0],
                        "lnum":inspect.getsourcelines(memb[1])[1],
                        "path":inspect.getsourcefile(memb[1])
                        })
          
    return {"astrotype":astrotype,
             "descriptorDict":ret}
