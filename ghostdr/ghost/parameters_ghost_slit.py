# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost_slit.py file, in alphabetical order.

from gempy.library import config
from geminidr.core import parameters_preprocess

class darkCorrectConfig(parameters_preprocess.darkCorrectConfig):
    def setDefaults(self):
        self.suffix = "_darkCorrected"
        self.dark = None
        self.do_cal = "skip"

class CRCorrectConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_CRCorrected",
                          optional=True)
    ndev = config.RangeField("Number of median absolute deviations for "
                             "CR identification", float, 20, min=3)
    max_iters = config.RangeField("Maximum number of iterations", int, 1, min=0)


class processSlitsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slitsProcessed",
                          optional=True)


class stackFramesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_stack", optional=True)
    create_multiple_stacks = config.Field(
        "Create a stack for each spectrograph observation? (ignored for calibrations)", bool, True)
