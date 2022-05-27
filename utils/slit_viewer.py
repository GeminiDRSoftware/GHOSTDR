#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import astropy.io.fits as pf
from matplotlib import pyplot as plt
from matplotlib import gridspec
from ghostdr.ghost.lookups import polyfit_dict
from ghostdr.ghost.polyfit import SlitView as sv
import astrodata
import ghost_instruments

# pylint: disable=invalid-name

parser = argparse.ArgumentParser()
parser.add_argument('path')
parser.add_argument('-save', action='store_true', help='save the '
                    'red and blue slit viewer cutouts into files '
                    'adjacent to the parent file')
parser.add_argument('-no_simult', action='store_true', help='remove '
                    'the simultaneous arc fibre, if any, before '
                    'generating the plots (only valid for High res '
                    'OBJECT frames, ignored otherwise)')
args = parser.parse_args()

try:
    hdus = pf.open(args.path)
except:
    print("Cannot open file: " + args.path)
    sys.exit()

try:
    res = 'high' if hdus[0].header['SMPNAME'] == 'HI_ONLY' else 'std'
except:
    print("File resolution cannot be determined")
    sys.exit()

try:
    objtype = hdus[0].header['OBSTYPE'] == 'OBJECT'
except:
    objtype = False

# Just accumulate the SCI extensions (and any unlabeled extensions)
scihdus = []
for h in hdus[1:]:
    try:
        if h.header['EXTNAME'] == 'SCI':
            scihdus.append(h)
    except KeyError:
        scihdus.append(h)

ad = astrodata.open(args.path)
res_mode = ad.res_mode()
binning=ad.detector_x_bin()
try:
    slitv_fn = polyfit_dict.get_polyfit_filename(None, 'slitv', res_mode, ad.ut_date(),
                                                 ad.filename, 'slitvmod')
    slitvpars = astrodata.open(slitv_fn)
    print(f"Using slitvmod {slitv_fn}")
except IOError:
    sys.exit(1)

# get the 2d profiles/images
red2d = [sv(h.data, None, slitvpars.TABLE[0], mode=res, binning=binning).cutout('red')
         for h in scihdus]
blu2d = [sv(h.data, None, slitvpars.TABLE[0], mode=res, binning=binning).cutout('blue')
         for h in scihdus]

# zero the simultaneous arc fibre pixels if requested and
# file is a high res object frame
if args.no_simult and res == 'high' and objtype:
    svo = sv(h.data, None, slitvpars.TABLE[0], mode=res, binning=binning)
    # get location of simultaneous arc (identically located
    # regardless of arm so just use red)
    bounds = svo.object_boundaries['red'][-1]
    for c in red2d+blu2d:
        c[bounds[0]:bounds[1]+1][:] = 0

# compute the 1d profiles
red1d = [np.sum(c, axis=1) for c in red2d]
blu1d = [np.sum(c, axis=1) for c in blu2d]

# build the composite image (of all slits side-by-side)
reds = red2d[0]
blues = blu2d[0]
for r, b in zip(red2d[1:], blu2d[1:]):
    reds = np.hstack((reds, r))
    blues = np.hstack((blues, b))

# generate the graphs/plots
for arm, im, lines in zip(['red', 'blue'], [reds, blues], [red1d, blu1d]):
    nlines = len(lines)
    fig = plt.figure(res+'/'+arm, figsize=(0.5+1.5*nlines, 6))
    gs = gridspec.GridSpec(1, 1+nlines, width_ratios=[0.5+1.2*nlines]+[2.5]*nlines)

    # image goes to the left
    axp = plt.subplot(gs[0])
    axp.imshow(im, origin='lower', aspect=1) #'auto')
    axp.set(ylim=[0, im.shape[0]])
    axp.get_xaxis().set_visible(False)

    # line profiles go to the right
    ax1 = plt.subplot(gs[1])
    ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
    axs = [plt.subplot(gs[i+2], sharex=ax1) for i in np.arange(nlines-1)]
    axs.insert(0, ax1)
    for ax in axs:
        labels = [item.get_text() for item in ax.get_yticklabels()]
        ax.set_yticklabels(['']*len(labels))
        ax.set(ylim=[0, im.shape[0]])
    [axs[i].plot(d, range(len(d)))  # pylint: disable=W0106
        for i, d in zip(np.arange(nlines), lines)]

    plt.subplots_adjust(wspace=0.3)
    fig.autofmt_xdate(rotation=75)
    gs.tight_layout(fig)

plt.show()

# save composite images (if requested) adjacent to the source file
if args.save:
    parts = os.path.splitext(args.path)
    for arm, im in zip(['_red', '_blue'], [reds, blues]):
        zlist = pf.HDUList(pf.PrimaryHDU())
        zlist.append(pf.ImageHDU(data=im))
        path = parts[0] + arm + parts[1]
        zlist.writeto(path, overwrite=True)
