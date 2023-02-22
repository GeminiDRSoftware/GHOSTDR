#   Copyright RSAA, The Australian National University
#   The Gemini High Resolution Optical Spectrograph
#   Project partners:
#       The Gemini Observatory (Gemini, AURA),
#       The Australian Astrophysical Observatory (AAO),
#       The Australian National University (ANU) and
#       The Canadian National Research Council (NRC)
#
# Script to calculate the fwhm from a raw slit viewer image.
# This script requires the sv_mask.fits calibration as generated from the
# process_sv_image.py script that is part of the instrument control system
# software.
#
# The fwhm calculated by this script is just an estimate.
# The input frame # is not bias/overscan/flat corrected at all.
#
# Usage:  python calculate_fwhm.py sv_image.fits
#   Args:
#    --binning - Binning of input image"  nargs=2, default=[1, 1], type=int
#    --mask - Slit viewer mask calibration", default="sv_mask.fits", type=str
#    fname - Input FITS filename  type=str
#
# Author: Jon Nielsen, ANU
#

import sys
import math
import re
import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from PIL import Image, ImageDraw

S32 = math.sqrt(3)/2.0
doprint = False
doplot = False

def getccdxy(ccd):
    m = re.search(r'(\d+):(\d+),(\d+):(\d+)', ccd)
    x = (int(m[1]) - 1, int(m[2]))
    y = (int(m[3]) - 1, int(m[4]))
    return (x, y)

class Fibers(object):
    """ Represent the details of a fiber bundle. """

    # The order in which fibers are presented at the slit
    fiber_order = []

    # The IFU that each fiber belongs to
    ifu = []

    # The focal plane offsets in x and y of each fiber
    offsets = {}

    um_per_arcsec = 610.0

    def __init__(self, lenslet_width):
        self.lenslet_width = lenslet_width

    def xyoffsets(self):
        """ Return the xy offsets of the fibers """
        #points = np.zeros((len(self.offsets), 3))
        points = {}
        for i,fiber in enumerate(self.offsets):
            #points[i] = self.offsets[fiber] + (fiber, )
            points[fiber] = 1000 * np.array(self.offsets[fiber]) * self.lenslet_width / self.um_per_arcsec
        # Convert to mas from number of fibers
        #points[:,0:2] = 1000 * (points[:,0:2] * self.lenslet_width/self.um_per_arcsec)
        return points

class SRFibers(Fibers):
    """ Represent the details of the standard resolution fibers.

    The fiber_order is the AAO's numbering system for the fibers in CY_RPT_50, where
    numbers in the first bundle (low and high res) go from 1 to 41, and the numbers in the
    second bundle go from 42 to 61. Fiber 62 is the simultaneous calibration fiber."""
    fiber_order = np.array([2, 5, 3, 1, 6, 4, 7, 14, 15, 16, 43, 46, 44, 42, 47, 45, 48])
    ifu = np.array([1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2])
    offsets = {
        # IFU 1 fibers
        1: (0, 0),
        2: (S32, 0.5),
        3: (S32, -0.5),
        4: (0, -1),
        5: (-S32, -0.5),
        6: (-S32, 0.5),
        7: (0, 1),
        # IFU 1 guide fibers
        8: (S32, 1.5),
        9: (2*S32, 0),
        10: (S32, -1.5),
        11: (-S32, -1.5),
        12: (-2*S32, 0),
        13: (-S32, 1.5),
        # Sky fibers 14, 15, 16
        14: (-S32/2, 0.5),
        15: (S32/2, 0),
        16: (-S32/2, -0.5),
        # IFU 2 fibers
        42: (0, 0),
        43: (S32, 0.5),
        44: (S32, -0.5),
        45: (0, -1),
        46: (-S32, -0.5),
        47: (-S32, 0.5),
        48: (0, 1),
        # IFU 2 guide fibers
        49: (S32, 1.5),
        50: (2*S32, 0),
        51: (S32, -1.5),
        52: (-S32, -1.5),
        53: (-2*S32, 0),
        54: (-S32, 1.5)
    }

    def __init__(self):
        lenslet_width = 240.0   # Lenslet flat-to-flat in microns
        super(SRFibers, self).__init__(lenslet_width)

class HRFibers(Fibers):
    """ Represent the details of the high resolution fibers.

    See SRFibers for detail.
    """
    fiber_order = np.array([62, 0, 25, 31, 27, 32, 26, 30, 18, 21, 19, 17, 22, 20, 23, 28, 34, 24, 29, 33, 35,
                            56, 59, 57, 55, 60, 58, 61])
    ifu = np.array([2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    offsets = {
        17: (0, 0),
        18: (S32, 0.5),
        19: (S32, -0.5),
        20: (0, -1),
        21: (-S32, -0.5),
        22: (-S32, 0.5),
        23: (0, 1),
        24: (S32, 1.5),
        25: (2*S32, 1),
        26: (2*S32, 0),
        27: (2*S32, -1),
        28: (S32, -1.5),
        29: (0, -2),
        30: (-S32, -1.5),
        31: (-2*S32, -1),
        32: (-2*S32, 0),
        33: (-2*S32, 1),
        34: (-S32, 1.5),
        35: (0, 2),
        # Guide fibers
        36: (0, 3),
        37: (3*S32, 1.5),
        38: (3*S32, -1.5),
        39: (0, -3),
        40: (-3*S32, -1.5),
        41: (-3*S32, 1.5),
        # Sky fibers 55-61
        55: (0, 0),
        56: (S32, 0.5),
        57: (S32, -0.5),
        58: (0, -1),
        59: (-S32, -0.5),
        60: (-S32, 0.5),
        61: (0, 1),
        # And fiber 62 is the simultaneous calibration source
        62: (0, 0)
    }

    def __init__(self):
        lenslet_width = 144.0  # Lenslet flat-to-flat in microns
        super(HRFibers, self).__init__(lenslet_width)

def main():
    # Parse the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--binning", help="Binning of input image", nargs=2, default=[None, None], type=int)
    parser.add_argument("--mask", help="Slit viewer mask calibration", default="sv_mask.fits", type=str)
    parser.add_argument("fname", help="Input FITS filename", type=str)
    args = parser.parse_args()

    # Open the input file and extract the image data
    # It could be in the PHU or extension 1
    # No provision is made for any other extensions at this time
    fitsfile = fits.open(args.fname)
    fitsdata = np.array(fitsfile[0].data, dtype=np.double)
    if fitsdata.size <= 1:
        fitsdata = np.array(fitsfile[1].data, dtype=np.double)
    if len(fitsfile) > 1:
        fitshdr = fitsfile[1].header
    else:
        fitshdr = fitsfile[0].header
    if doprint:
        print('Shape of input data is', fitsdata.shape)

    if fitshdr.get('CCDSUM'):
        biny = int(fitshdr['CCDSUM'].split()[0])
        binx = int(fitshdr['CCDSUM'].split()[1])
        if doprint:
            print('Binning set from FITS header:', biny, binx)
    else:
        biny = binx = 1
        print('Binning set to default:', biny, binx)

    if (biny < 1) or (binx < 1):
        print('Invalid binning', biny, binx)
        sys.exit(1)

    # Median of input data
    pixmed = np.median(fitsdata)
    # Report the median because it gives the user an idea of
    # which slits are illuminated or not
    print('Median pixel value of input frame =', pixmed)

    # FIXME Could subtract the median here if required
    # fitsdata = fitsdata - pixmed

    # Zoom the input data so it has the same size as a 1x1 binned image
    data = zoom(fitsdata, (biny, binx), order=0)

    # Create a full size image
    if fitshdr.get('CCDSIZE'):
        x, y = getccdxy(fitshdr['CCDSIZE'])
    else:
        x = (0, data.shape[1])
        y = (0, data.shape[0])
    fitsdata = np.zeros((y[1] - y[0], x[1] - x[0]))

    # Figure out where the pixels fit in this image
    if fitshdr.get('CCDSEC'):
        if doprint:
            print('CCDSEC', fitshdr['CCDSEC'])
        x, y = getccdxy(fitshdr['CCDSEC'])
    else:
        x = (0, data.shape[1])
        y = (0, data.shape[0])
    fitsdata[y[0]:y[1], x[0]:x[1]] = data

    if doplot:
        plt.imshow(fitsdata, origin='lower')
        plt.show()

    result = fitsdata

    # Load the calibration image
    # FIXME needs error handling
    npf = fits.open(args.mask)
    npimg = [npf[1].data, npf[2].data]

    # Figure out the fiber geometry
    sr = SRFibers()
    hr = HRFibers()

    # The fibers for each IFU
    ifus = {
             'std1': sr.fiber_order[sr.ifu == 1],
             'std2': sr.fiber_order[sr.ifu == 2],
             'hi': hr.fiber_order[hr.ifu == 1],
           }

    # The xy offsets (in milliarcseconds) for each fiber
    alloffsets = {**sr.xyoffsets(), **hr.xyoffsets()}

    if doprint:
        print(alloffsets)

    # Iterate over the red and blue slit images
    # (they have separate calibrations)
    for j, arm in enumerate(['red', 'blue']):

        # These are the fiber id numbers, which come from the fiber image
        labels = np.unique(npimg[j])

        # We count the flux in each fiber
        fluxes = np.zeros(np.max(labels) + 1, dtype=float)
        pixfluxes = np.zeros_like(fluxes, dtype=float)

        if doprint:
            print('labels', labels, labels.shape)
            print('result', result.shape)
            print('fluxes', fluxes.shape)

        for l in labels:
            if l < 1:
                # This is the background, so don't count it
                continue
            # These are the pixels that are in this fiber
            pixargs = npimg[j] == l
            # Count the total flux in this fiber
            fluxes[l] = np.sum(result[pixargs])
            pixfluxes[l] = fluxes[l] / np.count_nonzero(pixargs)

        # Now calculate the seeing for each IFU
        for ifu in ifus:
            # This is an array of [x, y, flux] for each fiber in this IFU
            ifu_flux = np.array([[alloffsets[i][0], alloffsets[i][1], fluxes[i]] for i in ifus[ifu]])
            avgperpix = np.average([pixfluxes[i] for i in ifus[ifu]])
            if doprint:
                print(ifu)
                print(ifus[ifu])
                print(ifu_flux)

            # This is the position of the brightest fiber in the array
            amax = np.argmax(ifu_flux[:, 2])

            # We take the x and y of the brightest fiber as an estimate of the centre of the star
            xmax = ifu_flux[amax][0]
            ymax = ifu_flux[amax][1]

            # The distance of each fiber from the brightest one
            fdist = np.array([math.sqrt((x - xmax)**2 + (y - ymax)**2) for (x, y) in zip(ifu_flux[:,0], ifu_flux[:,1])])

            # The log of the fluxes - we use the trick that the log of a gaussian looks like a quadratic
            log_fluxes = np.log(ifu_flux[:, 2])

            # And we also use the simplification that we are fitting a quadratic that has its maximum
            # at xmax, ymax. This means that we simply want y = a*x**2 + c, with no b*x term.
            # Thus the equation is linear in x**2 so we can use a simple linear fitting algorithm..

            # Fit an order 1 polynomial
            a, c = np.polyfit(fdist**2, log_fluxes, 1)

            # We want to know the value of x when the flux is half the maximum. Because we took a log
            # of the flux this means we are looking for y = c - log2, which simplifies down to
            # x = sqrt(-log2/a)
            halfmax_x_sq = -math.log(2) / a
            # Just check that we're not going to sqrt a negative number
            if (halfmax_x_sq < 0):
                halfmax_x = 0
            else:
                halfmax_x = math.sqrt(halfmax_x_sq)

            # The halfmax_x can be turned into a seeing value by doubling it (to get the full width)
            # And turning it into arcseconds
            fwhm = 2 * halfmax_x / 1000
            if doprint:
                print('halfmax_x', halfmax_x)
                print('halfmax =', a * halfmax_x**2 + c)
            print(arm, ifu, 'mean pixel value =', avgperpix, ', fwhm =', fwhm)

if __name__ == '__main__':
    main()
