recipe_tags = set(['GHOST', 'CAL', 'FLAT'])

from ghostdr.ghost.recipes.qa import recipes_FLAT


def recipeFlatCreateMaster(p):

    return recipes_FLAT.makeProcessedFlat(p)

_default = recipeFlatCreateMaster
