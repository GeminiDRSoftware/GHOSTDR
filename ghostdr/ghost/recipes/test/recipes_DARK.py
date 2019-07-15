recipe_tags = set(['GHOST', 'CAL', 'DARK'])

from ghostdr.ghost.recipes.qa import recipes_DARK


def recipeDarkCreateMaster(p):

    return recipes_DARK.makeProcessedDark(p)

_default = recipeDarkCreateMaster
