# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost_slit.py file, in alphabetical order.

from gempy.library import config
from geminidr.core import parameters_stack


class CRCorrectConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_CRCorrected",
                          optional=True)


class processSlitsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slitsProcessed",
                          optional=True)
    flat = config.Field("Slit flat field", str, None, optional=True)


class stackFramesConfig(parameters_stack.stackFramesConfig):
    def setDefaults(self):
        self.suffix = "_stack"
        self.nhigh = 1
        self.nlow = 1
        self.operation = "mean"
        self.reject_method = "none"
        self.apply_dq = True
        self.zero = False
        self.scale = False
        self.separate_ext = True
        self.statsec = None
