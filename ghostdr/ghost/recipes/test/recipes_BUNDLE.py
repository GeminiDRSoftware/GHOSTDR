recipe_tags = set(['GHOST', 'BUNDLE', 'RAW', 'UNPREPARED', ])

from ghostdr.ghost.recipes.qa.recipes_BUNDLE import makeProcessedBundle


def recipeBundleTest(p):
    return makeProcessedBundle(p)

_default = recipeBundleTest
