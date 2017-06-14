#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import astropy.io.fits as pf
from matplotlib import pyplot as plt
from astrodata_GHOST.polyfit import SlitView as sv

parser = argparse.ArgumentParser()
parser.add_argument('path')
parser.add_argument('-save', action='store_true', help='save the '
                    'red and blue slit viewer cutouts into files '
                    'adjacent to the parent file')
args = parser.parse_args()

try:
    hdus = pf.open(args.path)
except:
    print "Cannot open file: " + args.path
    sys.exit()

try:
    res = 'high' if hdus[0].header['SMPNAME'] == 'HI_ONLY' else 'std'
except:
    print "File resolution cannot be determined"
    sys.exit()

svobj = sv(hdus[1].data, None, mode=res)
reds = svobj.cutout('red')
blues = svobj.cutout('blue')

for ext in hdus[2:]:
    svobj = sv(ext.data, None, mode=res)
    reds = np.hstack((reds, svobj.cutout('red')))
    blues = np.hstack((blues, svobj.cutout('blue')))

plt.figure(res + '/red')
plt.imshow(reds, origin='lower')
plt.figure(res + '/blue')
plt.imshow(blues, origin='lower')
plt.show()

if args.save:
    parts = os.path.splitext(args.path)

    # save red cutouts
    redlist = pf.HDUList(pf.PrimaryHDU())
    redlist.append(pf.ImageHDU(data=reds))
    path = parts[0] + '_red' + parts[1]
    redlist.writeto(path, clobber=True)

    # save blue cutouts
    blulist = pf.HDUList(pf.PrimaryHDU())
    blulist.append(pf.ImageHDU(data=blues))
    path = parts[0] + '_blue' + parts[1]
    blulist.writeto(path, clobber=True)
