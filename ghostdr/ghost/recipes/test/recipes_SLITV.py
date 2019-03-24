recipe_tags = set(['GHOST', 'SLITV', ])

from ghostdr.ghost.recipes.sq.recipes_SLITV import makeProcessedSlit


def recipeSlitTest(p):
    return makeProcessedSlit(p)


def recipeRetrieveSlitArcTest(p):
    p.getProcessedArc()
    return


default = recipeSlitTest

