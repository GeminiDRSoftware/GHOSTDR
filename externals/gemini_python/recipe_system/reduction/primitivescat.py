#
#                                                                 gemini_python
#
#
#                                                               primitivescat.py
# ------------------------------------------------------------------------------
# $Id: primitivescat.py 5142 2015-02-17 21:39:45Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5142 $'[11:-2]
__version_date__ = '$Date: 2015-02-18 08:39:45 +1100 (Wed, 18 Feb 2015) $'[7:-2]
# ------------------------------------------------------------------------------

class PrimitivesCatalog(object):
    def __init__(self):
        self.catdict = {}
        
    def add_primitive_set(self, package, primsetEntry=None, primsetPath=None):
        pdict = {}
        self.catdict.update({primsetEntry : pdict})
        pdict.update({"package":package, "path":primsetPath})
        return
            
    def get_primcat_dict(self, primsetEntry):
        if primsetEntry in self.catdict:
            return self.catdict[primsetEntry]
        else:
            return None
