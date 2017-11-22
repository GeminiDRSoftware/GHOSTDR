"""
Recipes available to data with tags ['GHOST', 'SPECT', 'PARTNER_CAL'].
Default is "reduce". SQ is identical to QA recipe.
"""
recipe_tags = set(['GHOST', 'SPECT', 'PARTNER_CAL'])

from ..qa.recipes_SPECT_PARTNERCAL import reducePCal

default = reducePCal
