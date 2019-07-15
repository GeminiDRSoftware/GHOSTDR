recipe_tags = set(['GHOST', 'SLITV', 'BIAS', 'CAL', ])

from ghostdr.ghost.recipes.qa.recipes_BIAS_SLITV import makeProcessedBias


def recipeSlitBiasTest(p):
    return makeProcessedBias(p)


_default = recipeSlitBiasTest
