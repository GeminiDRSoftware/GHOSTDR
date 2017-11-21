"""
Recipes available to data with tags: GHOST, BUNDLE, NxN, UNPREPARED
Default is "makeProcessedBundle".
"""
recipe_tags = set(['GHOST', 'BUNDLE', 'NxN', 'UNPREPARED'])

def makeProcessedBundle(p):
    """
    This recipe processes GHOST observation bundles.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """
    p.splitBundle()
    return

default = makeProcessedBundle
