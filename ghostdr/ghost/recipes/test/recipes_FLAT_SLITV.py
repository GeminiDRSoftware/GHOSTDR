recipe_tags = set(['GHOST', 'SLITV', 'FLAT', 'CAL', ])

from ghostdr.ghost.recipes.sq.recipes_FLAT_SLITV import makeProcessedSlitFlat


def recipeSlitFlatTest(p):
    return makeProcessedSlitFlat(p)


def recipeRetrieveSlitDarkTest(p):
    p.getProcessedDark()
    return


default = recipeSlitFlatTest

