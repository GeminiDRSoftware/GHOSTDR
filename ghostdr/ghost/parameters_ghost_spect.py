# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost_spect.py file, in alphabetical order.

from .parameters_ghost import ParametersGHOST

class ParametersGHOSTSpect(ParametersGHOST):


    addWavelengthSolution = {
        "suffix": "_wavelengthAdded",
    }
    applyFlatBPM = {
        "suffix"            : "_flatBPMApplied",
        "flat"              : None,
        "flat_stream"       : None,
        "write_result"      : True,
    }
    barycentricCorrect = {
        "suffix": "_barycentricCorrected",
        "correction_factor": None,
    }
    clipSigmaBPM = {
#        "suffix"            : "_clipSigmaBPMApplied",
        "sigma"             : 3.0,
        "bpm_value"         : 1,
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
    rejectCosmicRays = {
        "suffix"            : "_cosmicRaysRejected",
        "subsampling"       : 2,
        "sigma_lim"         : 15.0,
        "f_lim"             : 5.0,
        "n_steps"           : 1,
    }
    responseCorrect = {
        "suffix"            : "_responseCorrected",
        "skip"              : False,
    }
    tileAmplifiers = {
        "suffix"            : "_ampsTiled",
    }
