"""
Recipes available to data with tags ``['GHOST', 'CAL', 'BIAS']``.
Default is ``makeProcessedBias``.
"""
recipe_tags = set(['GHOST', 'CAL', 'BIAS'])

from ..qa.recipes_BIAS import makeProcessedBias

_default = makeProcessedBias
