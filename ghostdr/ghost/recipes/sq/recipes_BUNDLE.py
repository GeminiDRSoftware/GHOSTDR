"""
Recipes available to data with tags ``['GHOST', `BUNDLE`]``.
Recipes imported from :any:`qa.recipes_BUNDLE`.
"""
recipe_tags = set(['GHOST', 'BUNDLE', 'RAW', 'UNPREPARED'])

from ..qa.recipes_BUNDLE import makeProcessedBundle

_default = makeProcessedBundle
