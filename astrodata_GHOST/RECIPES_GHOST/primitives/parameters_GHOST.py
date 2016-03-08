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
        }
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
}

