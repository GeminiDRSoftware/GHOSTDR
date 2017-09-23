"""
Recipes available to data with tags ['GHOST', 'CAL', 'ARC'].
Default is "makeProcessedArc". SQ is identical to QA recipe.
"""
recipe_tags = set(['GHOST', 'CAL', 'ARC'])

from ..qa.recipes_ARC import makeProcessedArc

default = makeProcessedArc
