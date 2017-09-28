
#                                                                  gemini_python
#
#                                                     primitives_ghost_bundle.py
# ------------------------------------------------------------------------------
from .primitives_ghost import GHOST
from .parameters_ghost_bundle import ParametersGHOSTBundle

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------
@parameter_override
class GHOSTBundle(GHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOSTBundle level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set(["GEMINI", "GHOST", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(GHOSTBundle, self).__init__(adinputs, **kwargs)
        self.parameters = ParametersGHOSTBundle

    # The splitting primitive should go here