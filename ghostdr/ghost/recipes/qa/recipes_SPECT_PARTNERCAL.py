"""
Recipes available to data with tags ['GHOST', 'SPECT', 'PARTNER_CAL'].
Default is "reduce".
"""
recipe_tags = set(['GHOST', 'SPECT', 'PARTNER_CAL'])

def reducePCal(p):
    """
    This recipe processes GHOST telluric standard data.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.darkCorrect()
    p.tileArrays()
    #p.rejectCosmicRays()
    p.applyFlatBPM() # Bitwise combine the flat BPM with the current BPM
                     # for the data. This is necessary because the flat isn't
                     # subtracted in the classical sense - rather, it's profile
                     # is subtracted from the object profile. Therefore, we
                     # must apply the BPM of the flat to the object file
                     # separately, before we extract its profile.
    p.extractProfile()
    p.flatCorrect() # Need to write our own, NOT USE GMOS - extract the flat profile,
                    # then simple division
    return

default = reduce
