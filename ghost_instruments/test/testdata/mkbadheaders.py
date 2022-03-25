import numpy as np
from astropy.io import fits

base = fits.open("bias_2_MEF_2x2_slit.fits")
base[0].header.pop('DATALAB', None)
base[0].header.pop('CCDSUM', None)
base[0].header['ARCBEFOR'] = True
base[1].header['CCDSUM'] = 'a b'
base[1].header['EXPTIME'] = -1
base.writeto("badheaders1.fits", overwrite=True)
base[0].header['ARCBEFOR'] = False
base[0].header['TILEARRY'] = True
base.writeto("badheaders2.fits", overwrite=True)
