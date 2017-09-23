# This parameter file contains the parameters related to the primitives located
# in the primitives_calibdb_ghost.py file, in alphabetical order.

from geminidr.core.parameters_calibdb import ParametersCalibDB

class ParametersCalibDBGHOST(ParametersCalibDB):

    storeProcessedSlit = {
        "suffix"            : "_slit",
    }

    storeProcessedSlitBias = {
        "suffix"            : "_slitbias",
    }

    storeProcessedSlitDark = {
        "suffix"            : "_slitdark",
    }

    storeProcessedSlitFlat = {
        "suffix"            : "_slitflat",
    }

    storeProcessedWavefit = {
        "suffix"            : "_wmodPolyfit",
    }

    storeProcessedXmod = {
        "suffix"            : "_xmodPolyfit",
    }