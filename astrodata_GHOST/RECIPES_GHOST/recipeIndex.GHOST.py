# As necessary, assign a specific recipe to GHOST AstroDataTypes
# Generic recipes from the astrodata_Gemini package should be available,
# for example makeProcessedBias.

localAstroTypeRecipeIndex = {
                            "GHOST_BIAS"   : ["makeProcessedBiasG"],
                            "GHOST_DARK"   : ["makeProcessedDarkG"],
                            "GHOST_ARC"    : ["makeProcessedArcG"],
                            "GHOST_FLAT"   : ["makeProcessedFlatG"],
                            "GHOST_OBJECT" : ["reduceG"],
                            }
