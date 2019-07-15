recipe_tags = set(['GHOST', 'SLITV', 'DARK', 'CAL', ])

from ghostdr.ghost.recipes.sq.recipes_DARK_SLITV import makeProcessedSlitDark


def recipeSlitDarkTest(p):
    return makeProcessedSlitDark(p)


def recipeRetrieveSlitBiasTest(p):
    p.getProcessedBias()
    return


_default = recipeSlitDarkTest
