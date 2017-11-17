"""
Recipes available to data with tags ['GHOST', 'SPECT'].
Default is "makeProcessedArc". SQ is identical to QA recipe.
"""
recipe_tags = set(['GHOST', 'SPECT'])

from ..qa.recipes_SPECT import reduce

default = reduce
