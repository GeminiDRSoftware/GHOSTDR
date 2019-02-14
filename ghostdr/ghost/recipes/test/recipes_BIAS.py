recipe_tags = set(['GHOST', 'CAL', 'BIAS'])

from ghostdr.ghost.recipes.qa import recipes_BIAS


def recipeBiasRemoveOverscan(p):
    p.prepare()
    p.overscanCorrect()
    # Note no addVAR or addDQ - want to leave the data and headers pristine
    return


def recipeBiasCreateMaster(p):

    # Copied SQ recipe - can't import directly (need to write after OS correct)
    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.writeOutputs()
    # p.tileArrays()
    p.addToList(purpose="forStack")
    p.getList(purpose="forStack")
    p.stackFrames(operation="median", apply_dq=True, reject_method='none')
    p.clipSigmaBPM(bpm_value=1, sigma=5.0)
    p.storeProcessedBias()
    return

    # p.prepare()
    # p.addDQ()
    # p.addVAR(read_noise=True)
    # p.overscanCorrect()
    # p.writeOutputs()
    # p.stackFrames()
    # p.clipSigmaBPM(bpm_value=1, sigma=5.0)
    # return

default = recipeBiasCreateMaster
