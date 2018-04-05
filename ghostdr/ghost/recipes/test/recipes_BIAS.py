recipe_tags = set(['GHOST', 'CAL', 'BIAS'])


def recipeBiasRemoveOverscan(p):
    p.prepare()
    p.overscanCorrect()
    # Note no addVAR or addDQ - want to leave the data and headers pristine
    return


def recipeBiasCreateMaster(p):
    p.prepare()
    p.overscanCorrect()
    p.writeOutputs()
    p.stackFrames()
    return

# default = recipeBiasRemoveOverscan
