#!/usr/bin/env python
import os
import astrodata
import ghost_instruments
from ghostdr.ghost.primitives_ghost import GHOST

all_files = [f for f in list(os.listdir('.')) if f.endswith('.fits')]
#if len(all_files) != 60:
#    raise IOError("Wrong number of files")

# Add DATALAB to SLIT files and tweak the others
# bias        1-00n
# arc (std)   2-001
# arc (high) 10-001
# flat (std)  3-00n
# flat (high) 4-00n
# dark        5-00n
# sky (std)   6-001
# sky (high)  7-001
# obj (std)   8-001 (0.5) 8-002 (1.0)
# obj (high)  9-001 (0.5) 9-002 (1.0)
for f in all_files:
    if 'SLIT' not in f:
        ad = astrodata.open(f)
        tags = ad.tags
        datalab = list(ad.data_label())
        if 'FLAT' in tags:
            datalab[14] = ('3' if 'std' in f else '4')
        elif f.startswith('sky'):
            datalab[14] = ('6' if 'std' in f else '7')
        elif f.startswith('obj'): 
            datalab[14] = ('8' if 'std' in f else '9')
            datalab[-1] = ('1' if '0.5' in f else '2')
        elif 'ARC' in f and 'HIGH' in f:
            datalab = datalab[:14] + ['1', '0'] + datalab[15:]
        ad.phu['DATALAB'] = ''.join(datalab)
        ad.write(clobber=True)

for f in all_files:
    if 'SLIT' in f:
        fname_start = f[:f.index('_SLIT')]
        data_label = astrodata.open([g for g in all_files
                                     if g.startswith(fname_start) and 'SLIT' not in g][0]).phu['DATALAB']
        ad = astrodata.open(f)
        ad.phu['DATALAB'] = data_label
        ad.write(clobber=True)

# Bias frames don't need calibrations so ditch them now
all_files = filter(lambda f: not f.startswith('bias'), all_files)

# Not the most efficient way, but the clearest (I think)
for f in all_files:
    ad = astrodata.open(f)
    tags = ad.tags
    p = GHOST([ad])
    arm = ad.arm()
    res = ad.res_mode()
    seeing = '0.5' if '0.5' in f else '1.0'

    if 'SLIT' in f:
        # All have a bias
        p.addCalibration(caltype='processed_bias', calfile='calibrations/processed_bias/bias_1_SLIT_bias.fits')
        if 'DARK' not in tags:
            p.addCalibration(caltype='processed_dark', calfile='calibrations/processed_dark/dark95_1_SLIT_dark.fits')
            if 'FLAT' not in tags:
                p.addCalibration(caltype='processed_slitflat', calfile='calibrations/processed_slitflat/flat95_{}_1_SLIT_slitflat.fits'.format(res))
    else:
        p.addCalibration(caltype='processed_bias', calfile='calibrations/processed_bias/bias_1_1x1_{}_bias.fits'.format(arm))
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
