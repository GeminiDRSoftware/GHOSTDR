# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost_slit.py file, in alphabetical order.

from .parameters_ghost import ParametersGHOST

class ParametersGHOSTSlit(ParametersGHOST):

    CRCorrect = {
        "suffix"            : "_CRCorrected",
    }
    processSlits = {
        "flat"              : None,
        "suffix"            : "_slitsProcessed",
    }
    stackFrames = {
        "suffix"            : "_stack",
        "nhigh"             : 1,
        "nlow"              : 1,
        "operation"         : "average",
        "reject_method"     : "none",
        "apply_dq"          : True,
        "zero"              : False,
        "scale"             : False,
        "separate_ext"      : True,
        "statsec"           : None,
    }
