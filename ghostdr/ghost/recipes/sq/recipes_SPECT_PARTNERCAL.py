"""
Recipes available to data with tags ``['GHOST', 'SPECT', 'PARTNER_CAL']``.
Default is ``reducePCal``, which is imported from
:any:`qa.recipes_SPECT_PARTNERCAL`.
"""
recipe_tags = set(['GHOST', 'SPECT', 'PARTNER_CAL'])

from ..qa.recipes_SPECT_PARTNERCAL import reducePCal

_default = reducePCal
