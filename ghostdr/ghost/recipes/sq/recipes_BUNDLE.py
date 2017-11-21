"""
Recipes available to data with tags ['GHOST', 'BUNDLE', 'NxN', UNPREPARED].
Default is "makeProcessedBundle". SQ is identical to QA recipe.
"""
recipe_tags = set(['GHOST', 'BUNDLE', 'NxN', 'UNPREPARED'])

from ..qa.recipes_BUNDLE import makeProcessedBundle

default = makeProcessedBundle
