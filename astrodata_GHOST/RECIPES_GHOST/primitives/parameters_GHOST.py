# Input parameters default settings for the GHOST primitives, whether they
# are defined in primitives_GHOST.py or inherited.
#
# This is a dictionary of dictionaries of dictionaries (... not that bad.)

{"myScienceStep" : {
    "threshold" : {
        "default"        : 3.0,
        "type"           : "int",
        "recipeOverride" : True,
        "userOverride"   : True,
        "uiLevel"        : "UIBASIC",
        },
    "nextparam" : {
        "default"        : "somestring",
        "type"           : "str",
        "recipeOverride" : True,
        "userOverride"   : True,
        "uiLevel"        : "UIBASIC",
        },
    },
 "anotherPrimitive" : {
    "combine" : {
        "default"        : "sum",
        "type"           : "str",
        "recipeOverride" : True,
        "userOverride"   : True,
        "uiLevel"        : "UIBASIC",
        },
    },
"standardizeHeaders":{
    "suffix":{
        "default"       : "_headersStandardized",
        "type"          : "str",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    },
"standardizeStructure":{
    "suffix":{
        "default"       : "_structureStandardized",
        "type"          : "str",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "attach_mdf":{
        "default"       : False,
        "type"          : "bool",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "mdf":{
        "default"       : None,
        # No default type defined, since the mdf parameter could be a string or
        # an AstroData object
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    },
 "validateData":{
    "suffix":{
        "default"       : "_dataValidated",
        "type"          : "str",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "repair":{
        "default"       : False,
        "type"          : "bool",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    },
"rejectCosmicRays": {
    "suffix":{
        "default"       : "_cosmicRaysRejected",
        "type"          : "str",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "subsampling":{
        "default"       : 2,
        "type"          : "int",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "sigma_lim":{
        "default"       : 3.0,
        "type"          : "float",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "f_lim":{
        "default"       : 3.0,
        "type"          : "float",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    "n_steps":{
        "default"       : 2,
        "type"          : "int",
        "recipeOverride": True,
        "userOverride"  : True,
        "uiLevel"       : "UIBASIC",
        },
    },
}

