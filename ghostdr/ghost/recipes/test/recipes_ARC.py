recipe_tags = set(['GHOST', 'CAL', 'FLAT'])

from ghostdr.ghost.recipes.qa import recipes_ARC


def recipeArcCreateMaster(p):

    return recipes_ARC.makeProcessedArc(p)

default = recipeArcCreateMaster
