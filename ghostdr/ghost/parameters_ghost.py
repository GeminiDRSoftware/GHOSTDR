# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost.py file, in alphabetical order.

from geminidr.core.parameters_ccd import ParametersCCD
from geminidr.gemini.parameters_gemini import ParametersGemini
from .parameters_calibdb_ghost import ParametersCalibDBGHOST

class ParametersGHOST(ParametersGemini, ParametersCCD, ParametersCalibDBGHOST):

    applyFlatBPM = {
        "suffix"            : "_flatBPMApplied",
        "flat"              : None,
        "flat_stream"       : None,
        "write_result"      : True,
    }
    clipSigmaBPM = {
#        "suffix"            : "_clipSigmaBPMApplied",
        "sigma"             : 3.0,
        "bpm_value"         : 1,
    }
    correctSlitCosmics = {
        "suffix"            : "_slitCosmicsCorrected",
    }
    extractProfile = {
        "suffix"            : "_extractedProfile",
        "slit"              : None,
        "slitflat"          : None,
        "write_result"      : True,
    }
    findApertures = {
        "slitflat"          : None,
        "suffix"            : "_findAper",
    }
    fitWavelength = {
        "suffix"            : "_fitWavl",
    }
    flatCorrect = {
        "suffix"            : "_flatCorrected",
        "flat"              : None,
        "slit"              : None,
        "slitflat"          : None,
        "write_result"      : True,
    }
    overscanCorrect = {
        "suffix"            : "_overscanCorrected",
        "niterate"          : 2,
        "high_reject"       : 3.0,
        "low_reject"        : 3.0,
        "function"          : "polynomial",
        "nbiascontam"       : 4,
        "order"             : 0,
    }
    processSlits = {
        "flat"              : None,
        "suffix"            : "_slitsProcessed",
    }
    rejectCosmicRays = {
        "suffix"            : "_cosmicRaysRejected",
        "subsampling"       : 2,
        "sigma_lim"         : 15.0,
        "f_lim"             : 5.0,
        "n_steps"           : 1,
    }
    stackSlitFrames = {
        "operation"         : "average",
        "reject_method"     : "none",
    }
    tileAmplifiers = {
        "suffix"            : "_ampsTiled",
    }
