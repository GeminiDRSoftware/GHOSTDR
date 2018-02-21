"""
Recipes available to data with tags ['GHOST', 'CAL', 'ARC'].
Default is "makeProcessedArc".
"""
recipe_tags = set(['GHOST', 'CAL', 'ARC'])

def makeProcessedArc(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input arc images into a single stacked arc image. This output
    processed arc is stored on disk using storeProcessedArc and has a name
    equal to the name of the first input arc image with "_arc.fits" appended.
    The wavelength solution is also stored.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    #p.tileArrays()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    # TODO? p.ADUToElectrons()
    p.darkCorrect()
    #p.rejectCosmicRays(
    # )
    p.tileArrays()
    p.extractProfile(sky_correct=False)
    p.fitWavelength()
    p.storeProcessedArc()
    return

default = makeProcessedArc
