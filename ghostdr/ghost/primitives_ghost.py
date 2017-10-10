#
#                                                                  gemini_python
#
#                                                            primitives_ghost.py
# ------------------------------------------------------------------------------
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core.primitives_ccd import CCD
from .primitives_calibdb_ghost import CalibDBGHOST

from .parameters_ghost import ParametersGHOST

from .lookups import timestamp_keywords as ghost_stamps

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------
@parameter_override
class GHOST(Gemini, CCD, CalibDBGHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOST level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set()  # Cannot be assigned as a class

    def __init__(self, adinputs, **kwargs):
        super(GHOST, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'ghostdr.ghost.lookups'
        self.parameters = ParametersGHOST
        # Add GHOST-specific timestamp keywords
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)
