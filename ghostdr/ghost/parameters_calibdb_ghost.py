# This parameter file contains the parameters related to the primitives located
# in the primitives_calibdb_ghost.py file, in alphabetical order.

from gempy.library import config


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
