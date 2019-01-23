recipe_tags = set(['GHOST', 'SLITV', 'ARC', 'CAL', ])

from ghostdr.ghost.recipes.sq.recipes_ARC_SLITV import makeProcessedSlitArc


def recipeSlitArcTest(p):
    return makeProcessedSlitArc(p)


def recipeRetrieveSlitFlatTest(p):
    p.getProcessedFlat()
    return


default = recipeSlitArcTest

