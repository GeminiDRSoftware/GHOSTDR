#!/usr/bin/env python
import os
import astrodata
import ghost_instruments
from ghostdr.ghost.primitives_ghost import GHOST
import pickle

all_files = [f for f in list(os.listdir('.')) if f.endswith('.fits')]

# don't process hotpixel planes, cosmic ray planes, check files, or bundle files (MEFs)
all_files = filter(lambda f: 'HP' not in f, all_files)
all_files = filter(lambda f: 'CR' not in f, all_files)
all_files = filter(lambda f: 'chk' not in f, all_files)
all_files = filter(lambda f: 'MEF' not in f, all_files)

# Bias frames don't need calibrations so ditch them now
all_files = filter(lambda f: not f.startswith('bias'), all_files)

# Not the most efficient way, but the clearest (I think)
for f in all_files:
    ad = astrodata.open(f)
    print(f)
    tags = ad.tags
    p = GHOST([ad])
    arm = ad.arm()
    res = ad.res_mode()
    # For the purposes of assigning the correct calibrations to different
    # binnings, find the binning section of the file name and save it
    try:
        x_index = f.index('x')
        binning = f[x_index-1:x_index+2]
    except: pass
    seeing = '0.5' if '0.5' in f else '1.0'

    if 'SLIT' in f:
        # All have a bias
        p.addCalibration(caltype='processed_bias', calfile='calibrations/processed_bias/bias_1_SLIT_bias.fits')
        if 'DARK' not in tags:
            p.addCalibration(caltype='processed_dark', calfile='calibrations/processed_dark/dark95_1_SLIT_dark.fits')
            if 'FLAT' not in tags:
                p.addCalibration(caltype='processed_slitflat', calfile='calibrations/processed_slitflat/flat95_{}_1_SLIT_slitflat.fits'.format(res))
    else:
        p.addCalibration(caltype='processed_bias', calfile='calibrations/processed_bias/bias_1_{}_{}_bias.fits'.format(binning, arm))
        if 'DARK' not in tags:
            p.addCalibration(caltype='processed_dark', calfile='calibrations/processed_dark/dark95_1_1x1_{}_dark.fits'.format(arm))
            if 'FLAT' in tags:
                p.addCalibration(caltype='processed_slitflat', calfile='calibrations/processed_slitflat/flat95_{}_1_SLIT_slitflat.fits'.format(res))
            else:
                if 'ARC' in tags:
                    p.addCalibration(caltype='processed_slit', calfile='calibrations/processed_slit/arc95_{}_SLIT_slit.fits'.format(res))
                else:
                    p.addCalibration(caltype='processed_slit', calfile='calibrations/processed_slit/obj95_{}_{}_SLIT_slit.fits'.format(seeing, res))
                p.addCalibration(caltype='processed_flat', calfile='calibrations/processed_flat/flat95_{}_1_1x1_{}_flat.fits'.format(res, arm))
                p.addCalibration(caltype='processed_slitflat', calfile='calibrations/processed_slitflat/flat95_{}_1_SLIT_slitflat.fits'.format(res))

