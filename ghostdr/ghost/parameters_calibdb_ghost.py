# This parameter file contains the parameters related to the primitives located
# in the primitives_calibdb_ghost.py file, in alphabetical order.

from gempy.library import config
from geminidr.core import parameters_calibdb


class storeCalibrationConfig(parameters_calibdb.storeCalibrationConfig):
    caltype = config.ChoiceField("Type of calibration", str,
                                 allowed={"processed_arc": "processed ARC",
                                          "processed_bias": "procsessed BIAS",
                                          "processed_dark": "processed DARK",
                                          "processed_flat": "processed FLAT",
                                          "processed_fringe": "processed fringe",
                                          "processed_slitflat": "processed slit-viewer flat",
                                          "processed_slit": "processed slit-viewer",
                                          "bpm": "bad pixel mask"},
                                 optional=False)


class getCalibrationConfig(parameters_calibdb.getCalibrationConfig):
    caltype = config.ChoiceField("Type of calibration", str,
                                 allowed={"processed_arc": "processed ARC",
                                          "processed_bias": "procsessed BIAS",
                                          "processed_dark": "processed DARK",
                                          "processed_flat": "processed FLAT",
                                          "processed_fringe": "processed fringe",
                                          "processed_slitflat": "processed slit-viewer flat",
                                          "processed_slit": "processed slit-viewer",
                                          "bpm": "bad pixel mask"},
                                 optional=False)

class getProcessedArcConfig(parameters_calibdb.getProcessedArcConfig):
    howmany = config.Field("How many arcs to return", int, 1, optional=True)



class storeProcessedSlitConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slit",
                          optional=True)


class storeProcessedSlitBiasConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slitbias",
                          optional=True)


class storeProcessedSlitDarkConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slitdark",
                          optional=True)


class storeProcessedSlitFlatConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slitflat",
                          optional=True)


class storeProcessedWavefitConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_wmodPolyfit",
                          optional=True)


class storeProcessedXmodConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_xmodPolyfit",
                          optional=True)

class getProcessedSlitFlatConfig(config.Config):
    refresh = config.Field("Refresh existing calibration associations?", bool,
                           True)

class getProcessedSlitConfig(config.Config):
    refresh = config.Field("Refresh existing calibration associations?", bool,
                           True)
