"""
Recipes available to data with tags ``['GHOST', 'CAL', 'ARC']``.
Default is ``makeProcessedArc``, which is imported from the QA module.
"""
recipe_tags = set(['GHOST', 'CAL', 'ARC'])

from ..qa.recipes_ARC import makeProcessedArc

default = makeProcessedArc
