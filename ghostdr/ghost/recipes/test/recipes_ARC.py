recipe_tags = set(['GHOST', 'CAL', 'ARC'])

from ghostdr.ghost.recipes.qa import recipes_ARC


def recipeArcCreateMaster(p):

    return recipes_ARC.makeProcessedArc(p)

_default = recipeArcCreateMaster
