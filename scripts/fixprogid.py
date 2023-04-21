from astropy.io import fits
from glob import glob
import sys

files = sorted([f for arg in sys.argv[1:] for f in glob(arg)])
for f in files:
    hlist = fits.open(f)
    print(f"Fixing {f}")
    hlist[0].header['GEMPRGID'] = "GS-2022B-Q-001"
    hlist.writeto(f, overwrite=True)

