"""
Recipes available to data with tags ['GHOST', 'SPECT'].
Default is "reduce".
"""
recipe_tags = set(['GHOST', 'SPECT'])

def reduce(p):
    """
    This recipe processes GHOST science data.

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
    #p.darkCorrect()
    p.tileArrays()
    #p.rejectCosmicRays()
    # TODO: Understand this! BPM combination is automatic in AD arithmetic
    p.applyFlatBPM() # Bitwise combine the flat BPM with the current BPM
                     # for the data
    p.extractProfile()
    p.flatCorrect() # Need to write our own, NOT USE GMOS - extract the flat profile,
                    # then simple division
    #p.responseCorrect() # canned standard star correction to the point of
                     # flatCorrect, plus a model flux
                     # User should be able to turn off this step
                     # Possible option for telluric correction (so
                     # option to override canned file used)
    #p.addWavelengthSolution() # should be able to accept multiple input
                           # arcs, e.g. from start and end of night,
                           # and interpolate in time
    #p.barycentricCorrect() # trivial - multiply wavelength scale
    # TODO: define reduce() and reduce_nostack() maybe?
    #p.interpolateAndCombine() # Should factor this step into separate recipe
    return

default = reduce
