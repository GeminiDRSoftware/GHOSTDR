# Input parameters default settings for the GHOST primitives, whether they
# are defined in primitives_GHOST.py or inherited.
#
# This is a dictionary of dictionaries of dictionaries (... not that bad.)

{
    "correctSlitCosmics": {
        "suffix": {
            "default":        "_slitCosmicsCorrected",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "extractProfile": {
      "suffix": {
            "default":        "_extractedProfile",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "findApertures": {
        "suffix": {
            "default":        "_findAper",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "fitWavelength": {
        "suffix": {
            "default":        "_fitWavl",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "processSlits": {
        "suffix": {
            "default":        "_slitsProcessed",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "rejectCosmicRays": {
        "suffix": {
            "default":        "_cosmicRaysRejected",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        "subsampling": {
            "default":        2,
            "type":           "int",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        "sigma_lim": {
            "default":        15.0,
            "type":           "float",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        "f_lim": {
            "default":        5.0,
            "type":           "float",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        "n_steps": {
            "default":        5,
            "type":           "int",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "standardizeHeaders": {
        "suffix": {
            "default":        "_headersStandardized",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "standardizeStructure": {
        "suffix": {
            "default":        "_structureStandardized",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedPolyfit": {
        "suffix": {
            "default":        "_xmodPolyfit",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedSlit": {
        "suffix": {
            "default":        "_slit",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedSlitBias": {
        "suffix": {
            "default":        "_slitBias",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedSlitDark": {
        "suffix": {
            "default":        "_slitDark",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedSlitFlat": {
        "suffix": {
            "default":        "_slitFlat",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "storeProcessedWavefit": {
        "suffix": {
            "default":        "_wmodPolyfit",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "tileAmplifiers": {
        "suffix": {
            "default":        "_tileAmplifiers",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    "validateData": {
        "suffix": {
            "default":        "_dataValidated",
            "type":           "str",
            "recipeOverride": True,
            "userOverride":   True,
            "uiLevel":        "UIBASIC",
            },
        },
    }
