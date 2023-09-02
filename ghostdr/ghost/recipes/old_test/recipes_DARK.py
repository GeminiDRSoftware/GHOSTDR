recipe_tags = set(['GHOST', 'CAL', 'DARK'])

from ghostdr.ghost.recipes.qa import recipes_DARK


def recipeDarkBiasCorrect(p):
    """
    This recipe performs the minimal steps required to old_test the bias
    subtraction.

    Parameters
    ----------
    p

    Returns
    -------

    """

    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    return

def recipeDarkCreateMaster(p):

    return recipes_DARK.makeProcessedDark(p)

_default = recipeDarkCreateMaster
