"""
Recipes available to data with tags ``['GHOST', 'SPECT']``.
Default is ``reduce``. SQ is identical to QA recipe, and is imported
from there.
"""
recipe_tags = set(['GHOST', 'SPECT'])

from ..qa.recipes_SPECT import reduce

default = reduce
