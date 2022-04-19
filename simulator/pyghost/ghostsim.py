"""This is a simple simulation code for GHOST or Veloce, with a class ARM that
simulates a single arm of the instrument. The key default parameters are
hardwired for each named configuration in the __init__ function of ARM.

Note that in this simulation code, the 'x' and 'y' directions are the
along-slit and dispersion directions respectively... (similar to physical
axes) but by convention, images are returned/displayed with a vertical slit
and a horizontal dispersion direction.

For a simple simulation, run:

import pyghost

blue = pyghost.ghostsim.Arm('blue')

blue.simulate_frame()
"""

from __future__ import print_function

import math
import os
# import pdb
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import numpy as np
import pylab as plt
import sys
import datetime
from dateutil import tz
from astropy import wcs
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import scipy.sparse
from copy import deepcopy

from pyghost import optics
from pyghost import cosmic

# uncomment to download the latest IERS predictions to avoid runtime warnings
from astropy.utils.data import download_file
from astropy.utils import iers
iers.IERS.iers_table = iers.IERS_A.open(
    download_file(iers.IERS_A_URL, cache=True))

try:
    import astropy.io.fits as pf
except ImportError:
    import pyfits as pf

# The directory that this file is located in
LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))

# Useful constants
PLANCK_H = 6.6256e-27  # Planck constant [erg s]
LIGHT_C = 2.99792458e18  # Speed of light [A/s]
S32 = math.sqrt(3)/2.0

# Number of fibers in each mode
# HR = high resolution
# SR = standard resolution
N_HR_SCI = 19
N_HR_SKY = 7
# We add 2, one for the ThAr fiber, and one for the gap
N_HR_TOT = N_HR_SCI + N_HR_SKY + 2
N_SR_SCI = 7
N_SR_SKY = 3
N_SR_TOT = 2 * N_SR_SCI + N_SR_SKY
N_GD = 6


def split_image(image, namps, return_headers=False):
    """ Split the input image into sections for each readout amplifier. """
    images = np.array_split(image, namps[0])
    current_size = 0
    ccdsecx = np.zeros((namps[1], 2))
    ccdsecy = np.zeros((namps[0], 2))
    newimages = []
    for i, im_amp in enumerate(images):
        ccdsecy[i, 0] = current_size
        current_size += im_amp.shape[0]
        ccdsecy[i, 1] = current_size
        newimages.extend(np.array_split(im_amp, namps[1], axis=1))
    images = newimages
    current_size = 0
    for i, im_amp in enumerate(images[0:namps[1]]):
        ccdsecx[i, 0] = current_size
        current_size += im_amp.shape[1]
        ccdsecx[i, 1] = current_size

    cxl, cyl = np.meshgrid(ccdsecx[:, 0], ccdsecy[:, 0])
    cxh, cyh = np.meshgrid(ccdsecx[:, 1], ccdsecy[:, 1])
    cxl = cxl.flatten()
    cyl = cyl.flatten()
    cxh = cxh.flatten()
    cyh = cyh.flatten()
    if return_headers:
        return images, cxl, cxh, cyl, cyh
    else:
        return images


def fftnoise(samples):
    """ Use an inverse FFT to generate noise. """
    samples = np.array(samples, dtype='complex')
    npoints = (len(samples) - 1) // 2
    phases = np.random.rand(npoints) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    samples[1:npoints+1] *= phases
    samples[-1:-1-npoints:-1] = np.conj(samples[1:npoints+1])
    return np.fft.ifft(samples).real


def frequency_noise(noise_freqs, sample_rate, shape, mean=0.0, std=1.0):
    """ Simulate noise at specific frequencies in a 1D array. """
    if np.isscalar(noise_freqs):
        noise_freqs = [noise_freqs]

    nsamples = np.prod(shape)
    sample_freqs = np.abs(np.fft.fftfreq(nsamples, 1/sample_rate))
    samples = np.zeros(nsamples)
    for freq in noise_freqs:
        idx = (np.abs(sample_freqs-freq)).argmin()
        samples[idx] = 1
    return mean + std/math.sqrt(2)*nsamples*fftnoise(samples).reshape(shape, order='F')


def thar_spectrum(ar_only=False):
    """Calculates a ThAr spectrum. Note that the flux scaling here is roughly correct
    for the lamp with no neutral density.

    Returns
    -------
    wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
    """

    if ar_only:
        thar = np.loadtxt(
            os.path.join(LOCAL_DIR, 'data/mnras_ar_only.txt'),
            usecols=[0, 1, 2])
    else:
        thar = np.loadtxt(
            os.path.join(LOCAL_DIR, 'data/mnras0378-0221-SD1.txt'),
            usecols=[0, 1, 2])
    # Create a fixed wavelength scale evenly spaced in log.
    thar_wave = 3600 * np.exp(np.arange(5e5)/5e5)
    thar_flux = np.zeros(int(5e5))
    # NB This is *not* perfect: we just find the nearest index of the
    # wavelength scale corresponding to each Th/Ar line.
    wave_ix = (np.log(thar[:, 1]/3600) * 5e5).astype(int)
    wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5-1).astype(int)
    thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))
    thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                            mode='same')
    # McDermid email of 10 July 2015 to Mike Ireland:
    # "For the arc, I integrated a single emission line of average brightness, and
    # putting this into a single resolution element in GHOST, get 370 e- per second per
    # pixel, averaging across the 3.4 pixel resolution element. "
    # This results in a peak flux of 370 * 3.4**2 * 7 = 30,000 for an "average" line
    # (of "strength" 3.0), or 300,000 for the brightest line (with a "strength" 4.0).
    thar_flux /= np.max(thar_flux) / 3e5
    return np.array([thar_wave/1e4, thar_flux])


def fits_in_dir(dirname, prefix=""):
    """ Return all fits files in the given directory """
    for fname in os.listdir(dirname):
        fullname = os.path.join(dirname, fname)
        if os.path.isfile(fullname) and fullname.endswith(".fits") and fname.startswith(prefix):
            yield fullname


def load_sky_from_dir(dirname):
    """ Load the UVES sky spectrum from the given directory """
    # Initialise our data structures
    wavel = []
    flux = []

    # Iterate over all the FITS files in the given directory
    for filename in fits_in_dir(dirname, prefix="fluxed_sky"):
        # Load the FITS hdulist
        hdulist = pf.open(filename)

        # Parse the WCS keywords in the primary HDU
        wav_wcs = wcs.WCS(hdulist[0].header)

        # Create an array of pixel coordinates
        pixcrd = np.array(np.arange(hdulist[0].data.shape[0]))

        # Convert pixel coordinates to world coordinates
        world = wav_wcs.wcs_pix2world(pixcrd, 0)[0]

        # Grab the data
        data = hdulist[0].data

        # Be careful - there's some dud data in there
        if filename.endswith("564U.fits"):
            args = (world > 5700) * (world < 5900)
            world = world[args]
            data = data[args]
        elif filename.endswith("800U.fits"):
            args = (world > 8530) * (world < 8630)
            world = world[args]
            data = data[args]

        # Accumulate the data
        wavel.extend(world)
        flux.extend(data)

    # Make sure the data is sorted by wavelength, in case we read
    # the files in some other order
    wavel = np.asarray(wavel)
    args = np.argsort(wavel)
    wavel = wavel[args]
    flux = np.asarray(flux)[args]
    # Why is there negative flux?
    flux[flux < 0] = 0

    return wavel, flux


def apply_binning(img, binmode):
    """ Bin the image (only works for power-of-2 binning factors) """
    for _ in range(int(math.log(binmode[0], 2))):
        img = img[:, ::2] + img[:, 1::2]
    for _ in range(int(math.log(binmode[1], 2))):
        img = img[::2, :] + img[1::2, :]
    # img /= float(binmode[0] * binmode[1])
    return img


def to_ushort(img):
    """ Convert to unsigned short, and deal with saturation """
    saturation = np.iinfo(np.uint16).max
    img[img < 0] = 0
    img[img > saturation] = saturation
    return np.asarray(img, dtype=np.uint16)


def add_overscan(img, oscan):
    """ add a 0-valued overscan region of width oscan to an image """
    return np.hstack((img, np.zeros((img.shape[0], oscan))))


class Fibers(object):
    """ Represent the details of a fiber bundle. """

    # The order in which fibers are presented at the slit
    fiber_order = []

    # The IFU that each fiber belongs to
    ifu = []

    # The focal plane offsets in x and y of each fiber
    offsets = {}

    hex_scale = 1.15

    def __init__(self, lenslet_width, microns_pix):
        self.lenslet_width = lenslet_width
        self.microns_pix = microns_pix

    def xyoffsets(self, ifu):
        """ Return the xy offsets of the fibers of the given IFU. """
        args = (self.ifu == ifu)
        points = np.zeros((np.count_nonzero(args), 2))
        for i, fiber in enumerate(self.fiber_order[args]):
            points[i] = self.offsets[fiber]
        points = (points * self.lenslet_width / self.microns_pix /
                  self.hex_scale).astype(int)
        return points[:, 0], points[:, 1]

    def plot(self):
        """ Plot the positions of the fibers for this IFU """
        points = np.array([(v[0], v[1], k, i)
                           for (k, v), i in zip(self.offsets.items(), self.ifu)
                           if v is not None])
        plt.plot(points[:, 0], points[:, 1], 'x')
        for x, y, label, ifu in points:
            plt.annotate(int(label), xy=(x, y), xytext=(0, (ifu-1)*10),
                         textcoords='offset points')
        plt.show()


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
        # Sky fibers
        14: None,
        15: None,
        16: None,
        # IFU 2 fibers
        42: (0, 0),
        43: (S32, 0.5),
        44: (S32, -0.5),
        45: (0, -1),
        46: (-S32, -0.5),
        47: (-S32, 0.5),
        48: (0, 1)
    }

    def __init__(self, lenslet_width, microns_pix):
        super(SRFibers, self).__init__(lenslet_width, microns_pix)


class HRFibers(Fibers):
    """ Represent the details of the high resolution fibers.

    See SRFibers for detail.
    """
    fiber_order = np.array([62, 0, 25, 31, 27, 32, 26, 30, 18, 21, 19, 17, 22, 20, 23, 28, 34, 24, 29, 33, 35,
                            56, 59, 57, 55, 60, 58, 61])
    ifu = np.array([0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
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
        # The following are the sky fibers
        55: None,
        56: None,
        57: None,
        58: None,
        59: None,
        60: None,
        61: None,
        # And fiber 62 is the simultaneous calibration source
        62: None
    }

    def __init__(self, lenslet_width, microns_pix):
        super(HRFibers, self).__init__(lenslet_width, microns_pix)


class SlitViewer(object):
    """A class for the slit viewer camera.  The initialisation function
    takes the duration of each slit viewer exposure, and the total number
    of exposures to be taken.  The data for each exposure will be filled
    in by each Arm object in the simulate_frame function.

    Parameters
    ----------
    cosmics: bool
        whether to add in cosmic rays or not

    crplane: bool
        Output a fits file containing locations where CRs were injected?

    split: bool
        Indicate whether to produce a single MEF containing the frames taken
        during the observation, or to generate individual files per frame
    """

    def __init__(self, cosmics, crplane, split):
        # FIXME these values are duplicated in the Arm class
        self.lenslet_std_size = 197.0   # Lenslet flat-to-flat in microns
        self.microns_pix = 2.0  # slit image microns per pixel
        # Below comes from CY_RPT_044, with 50mm/180mm demagnification and
        # the Bigeye G-283 (Sony ICX674) with 4.54 micron pixels.
        self.det_ysz = 1928  # TODO: per qsimaging.com, this is 1940
        self.det_xsz = 1452  # TODO: per qsimaging.com, this is 1460
        self.binning = 2  # Assumed that we bin in both directions
        # Effective pixel size of the slit viewer in microns.
        self.slitview_pxsize = 16.344 * self.binning
        self.slit_length = 3540.0    # Slit length in microns
        self.R_per_pixel = 195000.0  # Resolving power per pixel.
        self.slitcam_xsz = 160       # Slit viewer x size in binned pix.
        self.slitcam_ysz = 160       # Slit viewer y size in binned pix.
        self.slit_mode_offsets = {"high": -1000.0, "std": 1000.0}
        self.duration = 0
        self.sci_duration = None
        self.nexp = 1
        self.images = None
        self.cosims = []
        self.flux_profile = None
        self.utstart = datetime.datetime.utcnow()
        self.cosmics = cosmics
        self.crplane = crplane
        self.split = split

    def set_exposure(self, duration, sci_duration, flux_profile, utstart):
        """
        Set the exposure parameters for the slit viewing camera.

        Parameters
        ----------
        duration: float
            The duration of each slit viewing exposure

        nexp: int
            The number of slit viewing exposures to take

        flux_profile: int[2,n]
            The distribution of flux over the exposures.  Can be an array
            of any length - it will be interpolated and scaled to fit the
            actual number of exposures taken.

        utstart: datetime
            The UT start time from which to generate all time-related header
            keywords
        """
        self.duration = duration
        self.sci_duration = sci_duration
        if duration:  # to avoid divide by 0 (which uses nexp default of 1)
            self.nexp = int(math.ceil(float(sci_duration)/duration)) or 1
        self.utstart = utstart
        if flux_profile is None:
            self.flux_profile = np.ones((self.nexp))
        else:
            flux_profile = np.asarray(flux_profile)
            if (flux_profile[1].min() < 0) or (flux_profile[1].max() > 1):
                raise ValueError(
                    'Slit-viewing flux profile must have range in [0, 1]')
            # Ensure that the flux profile is evenly spaced over the number of
            # exposures we are taking.
            x = flux_profile[0] * self.nexp / flux_profile[0, -1]
            self.flux_profile = np.interp(
                np.arange(self.nexp), x, flux_profile[1])

        self.images = np.zeros(
            (self.nexp, self.slitcam_ysz, self.slitcam_xsz), dtype=int)
        if self.cosmics:
            self.cosims = [
                cosmic.cosmic(tuple(self.binning*x for x in im.shape), duration,
                10, 2.0, False, [4.54, 4.54, 10]) for im in self.images]

    def save(self, fname, obstype, res, bias=1000, readnoise=8.0, gain=1.0,
             data_label=1):
        """
        Save the slit viewer frames in self.images.

        Parameters
        ----------
        fname: string
            FITS filename for saving the images
        """

        if self.images is None:
            return

        # Calculate the detector/ccd section
        y0 = self.det_ysz // 2 - self.slitcam_ysz
        x0 = self.det_xsz // 2 - self.slitcam_xsz
        secstr = '[{}:{},{}:{}]'.format(
            y0, y0 + self.slitcam_ysz*self.binning - 1,
            x0, x0 + self.slitcam_xsz*self.binning - 1)
        detsz = '[1:{},1:{}]'.format(self.det_ysz, self.det_xsz)

        header = pf.Header()
        header['CAMERA'] = ('Slit', 'Camera name')
        header['OBSERVAT'] = (
            'Gemini-South', 'Name of telescope (Gemini-North|Gemini-South)')
        header['TELESCOP'] = 'Gemini-South'
        header['INSTRUME'] = ('GHOST', 'Instrument used to acquire data')
        if obstype=='STANDARD':
            header['OBSTYPE'] = ('OBJECT', 'Observation type')
        else:
            header['OBSTYPE'] = (obstype, 'Observation type')
        header['CCDNAME'] = ('Sony-ICX674', 'CCD name')
        header['ORIGFN'] = fname + 'SLIT.fits'

        # required by the calibration manager
        header['RAWPIREQ'] = 'yes'  # 'no', 'check', 'unknown'
        header['RAWGEMQA'] = 'usable'  # 'bad', 'check', 'unknown'

        # populate OBSCLASS keyword
        obsclass = dict(FLAT='partnerCal', ARC='partnerCal',  # noqa
                        OBJECT='science', BIAS='dayCal', DARK='dayCal', SKY='',
                        STANDARD='partnerCal')
        header['OBSCLASS'] = (obsclass[obstype], 'Observe class')

        # populate GEMPRGID, OBSID, and DATALAB keywords
        prgid = 'GS-2016B-Q-20'  # just choose a base at random
        header['GEMPRGID'] = (prgid, 'Gemini programme ID')
        obsid = dict(
            high=dict(BIAS='1', DARK='5', FLAT='4', ARC='10', OBJECT='9',
                      SKY='7', STANDARD='11'),  # noqa
            std=dict(BIAS='1', DARK='5', FLAT='3', ARC='2', OBJECT='8',
                     SKY='6', STANDARD='10'),  # noqa
            )
        header['OBSID'] = (
            prgid+'-'+obsid[res][obstype], 'Observation ID / Data label')
        header['DATALAB'] = (prgid+'-'+obsid[res][obstype]+'-' + (
            '%03d' % data_label), 'DHS data label')

        # perhaps wrongly (and out of line with their comments), the ghost
        # recipes interpret the following as the start/end times of the science
        # exposure taken concurrently with this series of slit viewer images
        # (and when the MEF-splitter primitive is written, it must remember to
        # populate them as such)
        header['DATE-OBS'] = (
            self.utstart.strftime("%Y-%m-%d"), 'UT date at observation start')
        header['UTSTART'] = (self.utstart.strftime("%H:%M:%S.%f")[:-3],
            'UT time at observation start')  # noqa
        header['UTEND'] = (
            (self.utstart + datetime.timedelta(seconds=self.sci_duration))
            .strftime("%H:%M:%S.%f")[:-3], 'UT time at observation end')

        # resolution-related keywords
        if obstype != 'BIAS' and obstype != 'DARK':
            header['SMPNAME'] = ('HI_ONLY' if res == 'high' else 'LO_ONLY')
            header['SMPPOS'] = (1 if res == 'high' else 2)

        hdulist = pf.HDUList([pf.PrimaryHDU(header=header)])
        crhdu = pf.HDUList(pf.PrimaryHDU(header=pf.Header()))
        for expid, image in enumerate(self.images):
            # Construct a new header
            hdr = pf.Header()

            cosim = self.cosims[expid]
            cosim = apply_binning(cosim, (self.binning, self.binning))
            image += to_ushort(cosim)

            # Grab the image data
            data = image / gain

            # Add some readout noise (linear regime assumed, so order doesn't
            # matter)
            data += np.random.normal(
                loc=bias, scale=readnoise / gain, size=data.shape)

            # Now restrict to 14 bits (BigEye G283) and digitize
            data = np.minimum(np.maximum(data, 0), 1 << 14).astype(np.int16)

            # Set the id of this extension
            hdr['EXPID'] = expid+1

            # And all the other headers
            hdr['EXTNAME'] = ('SCI', 'extension name')
            hdr['RDNOISE'] = (readnoise, 'Readout noise')
            hdr['GAIN'] = (gain, 'Amplifier gain')
            hdr['CAMERA'] = ('Slit', 'Camera name')
            hdr['DETSIZE'] = detsz
            hdr['CCDSIZE'] = detsz
            hdr['DETSEC'] = secstr
            hdr['CCDSEC'] = secstr
            datasec = '[1:{1},1:{0}]'.format(*data.shape)
            hdr['DATASEC'] = (datasec, 'Data section(s)')
            hdr['TRIMSEC'] = (datasec, 'Trim section(s)')
            hdr['AMPNAME'] = (0, 'Amplifier name')
            hdr['CCDNAME'] = ('Sony-ICX674', 'CCD name')
            hdr['CCDSUM'] = str(self.binning) + " " + str(self.binning)

            utcnow = self.utstart + datetime.timedelta(
                seconds=expid*self.duration)
            hdr['EXPUTST'] = (utcnow.strftime("%H:%M:%S.%f")[:-3],
                'UT time at exposure start')  # noqa
            hdr['EXPUTEND'] = (
                (utcnow + datetime.timedelta(seconds=self.duration))
                .strftime("%H:%M:%S.%f")[:-3], 'UT time at exposure end')
            # needed by gemcombine.cl (used by stackFrames())
            hdr['EXPTIME'] = self.duration

            if self.cosmics:
                cosim = self.cosims[expid]
                hdr['NCRPIX'] = np.count_nonzero(cosim)

                if self.crplane:
                    crhdr = pf.Header()
                    crhdr['DETSIZE'] = (detsz, 'Detector size')
                    crhdr['DETSEC'] = (secstr, 'Detector section(s)')
                    crhdu.append(
                        pf.ImageHDU(data=to_ushort(cosim), header=crhdr))

            hdulist.append(pf.ImageHDU(data=data, header=hdr))

        if not self.split:
            return hdulist

        print('Writing ' + fname + 'SLIT.fits')
        hdulist.writeto(fname + 'SLIT.fits', overwrite=True)
        if self.cosmics and self.crplane:
            print('Writing ' + fname + 'SLIT_CR.fits')
            crhdu.writeto(fname + 'SLIT_CR.fits', overwrite=True)

    def create_slitview_frames(self, im_slit, spectrum, slitview_wave,
                               slitview_frac, slitview_offset, mode='high'):
        """
        Create slit viewer frames, and store them in self.images.

        Parameters
        ----------
        im_slit: numpy array
            Image of the slit, sampled at self.microns_pix sampling

        spectrum: numpy array
            (2, nwave) array of wavelengths and fluxes

        slitview_wave: array (2)
            the low and high wavelength cutoffs of this exposure

        slitview_frac
            the fraction of flux that is diverted into the slit viewing camera

        slitview_offset
            the slit plane offset in microns

        mode: string
            'high' or 'std'

        """

        # If the images have not been set up then there are no exposures on this
        # camera.
        if self.images is None:
            return

        # Create a sub-array with just enough pixels to sample the slit.
        n_xpix = int(self.slit_length / self.slitview_pxsize) + 2
        n_ypix = int(self.lenslet_std_size / self.slitview_pxsize) + 2
        xy = np.meshgrid(np.arange(n_xpix) - n_xpix // 2,
                         np.arange(n_ypix) - n_ypix // 2)
        ycoord = (xy[1] * self.slitview_pxsize / self.microns_pix).astype(int)
        ycoord += im_slit.shape[0]//2
        xcoord = (xy[0] * self.slitview_pxsize / self.microns_pix).astype(int)
        xcoord += im_slit.shape[1]//2

        slit_camera_image = np.zeros((n_ypix, n_xpix))
        # Roughly convolve with a pixel with a 4-point dither.
        dither_pix = int(0.25 * self.slitview_pxsize / self.microns_pix)
        for xdither in [-dither_pix, dither_pix]:
            for ydither in [-dither_pix, dither_pix]:
                slit_camera_image += im_slit[ycoord + ydither, xcoord + xdither]

        # Re-normalize
        slit_camera_image /= 4
        slit_camera_image *= (self.slitview_pxsize / self.microns_pix)**2

        # Multiply by total flux within our bandpass. First multiply by mean
        # flux # per resolution element, then scale WARNING: Doesn't work if
        # min and max outside
        wave_in_filter = (slitview_wave[0] < spectrum[0]) & \
                         (spectrum[0] < slitview_wave[1])
        if sum(wave_in_filter) == 0:
            wave_in_filter = np.argmin(np.abs(
                spectrum[0] - 0.5 * slitview_wave[0] - 0.5 * slitview_wave[1]))
        slit_camera_image *= np.mean(spectrum[1, wave_in_filter])
        frac_bw = 0.5*(slitview_wave[1] - slitview_wave[0]) / \
            (slitview_wave[1] + slitview_wave[0])
        slit_camera_image *= self.duration * self.R_per_pixel * frac_bw * \
            slitview_frac
        slit_camera_image = np.maximum(slit_camera_image, 0)

        # TODO Determine if flipping is required once we have real data
        # MCW 190911 - Flip the image y-axis (axis 0)
        # slit_camera_image = slit_camera_image[::-1, :]

        # Store the full-sized image in units of photons.
        xoffset = int(slitview_offset[0]/self.slitview_pxsize)
        yoffset = int((self.slit_mode_offsets[mode] +
                       slitview_offset[1]) / self.slitview_pxsize)
        xy = np.meshgrid(
            np.arange(n_xpix) - n_xpix // 2 + xoffset + self.slitcam_xsz // 2,
            np.arange(n_ypix) - n_ypix // 2 + yoffset + self.slitcam_ysz // 2)

        for image, profile in zip(self.images, self.flux_profile):
            image[tuple(xy)] += np.random.poisson(profile * slit_camera_image)


class Arm(object):
    """A class for each arm of the spectrograph.

    Parameters
    ----------
    arm: string
        identifies which arm of the GHOST spectrograph. one of either 'red' or
        'blue'

    crplane: bool
        Output a fits file containing locations where CRs were injected?

    split: bool
        Indicate whether to produce a single MEF containing the frames taken
        during the observation, or to generate individual files per frame

    check: bool
        Output a fits file containing input spectra interpolated onto
        the pixel grid for each order?

    dual_target: bool
        Output STD res fits files that contain object spectra in both IFUs
        (Ignored for HIGH res mode)
    """

    ARM_OPTIONS = [
        'red',
        'blue',
    ]
    RES_OPTIONS = [
        'high',
        'std',
    ]

    # pylint: disable=too-many-instance-attributes

    def __init__(self, arm, slit_viewer=None, cosmics=True, crplane=False,
                 hpplane=False, split=False, check=False, dual_target=False):
        if arm.lower() not in self.ARM_OPTIONS:
            raise ValueError('arm must be one of %s' % (
                ','.join(self.ARM_OPTIONS),
            ))
        self.arm = arm.lower()
        self.d_y = 1000./52.67          # Distance in microns
        self.theta = 65.0               # Blaze angle
        self.assym = 1.0/0.41           # Magnification
        self.gamma = 0.56      # Echelle gamma
        self.nwave = int(1e2)       # Wavelengths per order for interpolation.
        self.f_col = 1750.6    # Collimator focal length.
        self.lenslet_high_size = 118.0  # Lenslet flat-to-flat in microns
        self.lenslet_std_size = 197.0   # Lenslet flat-to-flat in microns
        self.microns_pix = 2.0  # slit image microns per pixel
        self.microns_arcsec = 400.0  # slit image plane microns per arcsec
        self.im_slit_sz = 2048  # Size of the image slit size in pixels.
        self.sample_rate = 1e6  # Sample rate of the pixels
        self.cosmics = cosmics
        self.crplane = crplane
        self.hpplane = hpplane
        self.split = split
        self.check = check
        self.dual_target=dual_target
        self.slitv = slit_viewer
        # Do we need to re-compute the spectral format?
        self.stale_spectral_format = True
        # self.set_mode(mode)
        if self.arm == 'red':
            # Additional slit rotation across an order needed to match Zemax.
            self.extra_rot = 2.0
            self.szx = 6160        # Switched by MCWMCW
            self.szy = 6144        # 190814
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.drot = -2.0       # Detector rotation
            self.d_x = 1000/565.   # VPH line spacing
            self.theta_i = 30.0    # Prism incidence angle
            self.alpha1 = 0.0      # First prism apex angle
            self.alpha2 = 0.0      # Second prism apex angle
            self.order_min = 34
            self.order_max = 67
            self.dettype = 'E2V-CCD-231-C6'
            self.thick = 40        # detector thickness in microns
            self.slitview_wave = [0.6, 0.76]  # wavelength in microns
            self.slitview_frac = 0.02         # Fractional flux.
            self.slitview_offset = [15.0, -1500.0]  # (x,y) plane offset, micron
        elif self.arm == 'blue':
            # Additional slit rotation accross an order needed to match Zemax.
            self.extra_rot = 2.0
            self.szx = 4112  # Switched by MCW
            self.szy = 4096  # 190814
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.d_x = 1000/1137.  # VPH line spacing
            self.theta_i = 30.0    # Prism incidence angle
            self.drot = -2.0       # Detector rotation.
            self.alpha1 = 0.0      # First prism apex angle
            self.alpha2 = 0.0      # Second prism apex angle
            self.order_min = 63
            self.order_max = 95
            self.dettype = 'E2V-CCD-231-84'
            self.thick = 16        # detector thickness in microns
            self.slitview_wave = [0.43, 0.6]  # wavelength in microns
            self.slitview_frac = 0.02         # Fractional flux.
            self.slitview_offset = [-15.0, 1500.0]  # (x,y) plane offset, micron
        else:
            raise RuntimeError('Order information not provided in Arm class '
                               'for arm %s - aborting' % (self.arm, ))

        # define hot pixels for this detector (applied systematically to all
        # frames generated by the current run of this simulator)
        shape = (self.szy, self.szx)
        nbads = [3800, 5000][self.arm == 'red']
        idx = tuple([np.random.choice(dm, nbads) for dm in shape])
        vals = np.random.randint(200, 50000, nbads) * 1.0
        self.hot_pix = scipy.sparse.coo_matrix((vals, idx), shape=shape)

    def set_mode(self, new_mode):
        """Set a new mode (high or standard res) for tramline purposes.

        THIS IS CURRENTLY ONLY USED FOR EXTRACTION (PYMFE) TESTING

        A SIMILAR ROUTINE NEEDS TO START WITH HEADER KEYWORDS THEN THE
        FLUXES FROM THE SLIT VIEWING CAMERA.
        """
        self.mode = new_mode.lower()
        if self.mode == 'high':
            self.lenslet_width = self.lenslet_high_size
            self.nl = 28
            # Set default profiles - object, sky and reference
            fluxes = np.zeros((self.nl, 3))
            fluxes[2:21, 0] = 0.37
            fluxes[8:15, 0] = 0.78
            fluxes[11, 0] = 1.0
            # NB if on the following line, fluxes[2:,1]=1.0 is set, sky will be
            # subtracted automatically.
            fluxes[2+19:, 1] = 1.0
            fluxes[0, 2] = 1.0
        elif self.mode == 'std':
            self.lenslet_width = self.lenslet_std_size
            self.nl = 17
            # Set default profiles - object 1, sky and object 2
            fluxes = np.zeros((self.nl, 3))
            fluxes[0:7, 0]  = 1.0
            fluxes[7:10, 1] = 1.0
            fluxes[10:, 2] = 1.0
        else:
            print("Unknown mode!")
            raise UserWarning
        self.fluxes = fluxes

    def spectral_format(self, xoff=0.0, yoff=0.0, ccd_centre=None, verbose=False):
        """Create a spectrum, with wavelengths sampled in 2 orders.

        Parameters
        ----------
        xoff: float
            An input offset from the field center in the slit plane in
            mm in the x (spatial) direction.
        yoff: float
            An input offset from the field center in the slit plane in
            mm in the y (spectral) direction.
        ccd_centre: dict
            An input describing internal parameters for the angle of the
            center of the CCD. To run this program multiple times with the
            same co-ordinate system, take the returned ccd_centre and use it
            as an input.

        Returns
        -------
        x:  (norders, n_y) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel
            and each order (m).
        wave: (norders, n_y) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (norders, n_y) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        ccd_centre: dict
            Parameters of the internal co-ordinate system describing the center
            of the CCD.
        """
        # Parameters for the Echelle. Note that we put the
        # co-ordinate system along the principle Echelle axis, and
        # make the beam come in at the gamma angle.
        u_1 = -np.sin(np.radians(self.gamma) + xoff/self.f_col)
        u_2 = np.sin(yoff/self.f_col)
        u_3 = np.sqrt(1 - u_1**2 - u_2**2)
        u_vect = np.array([u_1, u_2, u_3])
        l_vect = np.array([1.0, 0, 0])
        s_vect = np.array([0, np.cos(np.radians(self.theta)),
                           -np.sin(np.radians(self.theta))])
        # Orders for each wavelength. We choose +/- 1 free spectral range.
        orders = np.arange(self.order_min, self.order_max+1, dtype=int)
        wave_mins = 2*self.d_y*np.sin(np.radians(self.theta))/(orders + 1.0)
        wave_maxs = 2*self.d_y*np.sin(np.radians(self.theta))/(orders - 1.0)
        wave = np.empty((len(orders), self.nwave))
        for i in range(len(orders)):
            wave[i, :] = np.linspace(wave_mins[i], wave_maxs[i], self.nwave)
        wave = wave.flatten()
        orders = np.repeat(orders, self.nwave)
        order_frac = np.abs(orders -
                            2*self.d_y*np.sin(np.radians(self.theta))/wave)
        ml_d = orders*wave/self.d_y
        # Propagate the beam through the Echelle.
        v_vects = np.zeros((3, len(wave)))
        for i in range(len(wave)):
            v_vects[:, i] = optics.grating_sim(u_vect, l_vect, s_vect, ml_d[i])
        # Find the current mean direction in the x-z plane, and magnify
        # the angles to represent passage through the beam reducer.
        if ccd_centre:
            mean_v = ccd_centre['mean_v']
        else:
            mean_v = np.mean(v_vects, axis=1)
            # As the range of angles is so large in the y direction, the mean
            # will depend on the wavelength sampling within an order. So just
            # consider a horizontal beam.
            mean_v[1] = 0
            # Re-normalise this mean direction vector
            mean_v /= np.sqrt(np.sum(mean_v**2))

        for i in range(len(wave)):
            # Expand the range of angles around the mean direction.
            temp = mean_v + (v_vects[:, i]-mean_v)*self.assym
            # Re-normalise.
            v_vects[:, i] = temp/np.sum(temp**2)

        # Here we diverge from Veloce. We will ignore the glass, and
        # just consider the cross-disperser.
        l_vect = np.array([0, -1, 0])
        theta_xdp = -self.theta_i + self.gamma
        # Angle on next line may be negative...
        s_vect = optics.rotate_xz(np.array([1, 0, 0]), theta_xdp)
        n_vect = np.cross(s_vect, l_vect)  # The normal
        incidence_angle = np.degrees(np.arccos(np.dot(mean_v, n_vect)))
        if verbose:
            print('Incidence angle in air: {0:5.3f}'.format(incidence_angle))
        # w is the exit vector after the grating.
        w_vects = np.zeros((3, len(wave)))
        for i in range(len(wave)):
            w_vects[:, i] = optics.grating_sim(v_vects[:, i], l_vect, s_vect,
                                               wave[i]/self.d_x)
        mean_w = np.mean(w_vects, axis=1)
        mean_w[1] = 0
        mean_w /= np.sqrt(np.sum(mean_w**2))
        exit_angle = np.degrees(np.arccos(np.dot(mean_w, n_vect)))
        if verbose:
            print('Grating exit angle in glass: {0:5.3f}'.format(exit_angle))
        # Define the CCD x and y axes by the spread of angles.
        if ccd_centre:
            ccdx = ccd_centre['ccdx']
            ccdy = ccd_centre['ccdy']
        else:
            ccdy = np.array([0, 1, 0])
            ccdx = np.array([1, 0, 0]) - np.dot([1, 0, 0], mean_w)*mean_w
            ccdx[1] = 0
            ccdx /= np.sqrt(np.sum(ccdx**2))

        # Make the spectrum on the detector.
        xpx = np.zeros(len(wave))
        ypx = np.zeros(len(wave))
        # There is definitely a more vectorised way to do this.
        for i in range(len(wave)):
            xpx[i] = np.dot(ccdx, w_vects[:, i])*self.f_cam/self.px_sz
            ypx[i] = np.dot(ccdy, w_vects[:, i])*self.f_cam/self.px_sz
            # Rotate the chip to get the orders along the columns.
            rot_rad = np.radians(self.drot)
            rot_matrix = np.array([[np.cos(rot_rad), np.sin(rot_rad)],
                                   [-np.sin(rot_rad), np.cos(rot_rad)]])
            [xpx[i], ypx[i]] = np.dot(rot_matrix, [xpx[i], ypx[i]])
        # Center the spectra on the CCD in the x-direction.
        if ccd_centre:
            xpix_offset = ccd_centre['xpix_offset']
        else:
            w_ix = np.where((ypx < self.szy//2) * (ypx > -self.szy//2))[0]
            xpix_offset = 0.5*(np.min(xpx[w_ix]) + np.max(xpx[w_ix]))

        xpx -= xpix_offset
        # Now lets interpolate onto a pixel grid rather than the
        # arbitrary wavelength grid we began with.
        n_orders = self.order_max-self.order_min+1
        x_int = np.zeros((n_orders, self.szy))
        wave_int = np.zeros((n_orders, self.szy))
        blaze_int = np.zeros((n_orders, self.szy))
        # plt.clf()
        for order in range(self.order_min, self.order_max+1):
            w_ix = np.where(orders == order)[0]
            ypx_min = np.max([np.min(ypx[w_ix]).astype(int), -self.szy//2])
            ypx_max = np.min([np.max(ypx[w_ix]).astype(int), self.szy//2])
            y_int_m = np.arange(ypx_min, ypx_max, dtype=int)
            y_ix = y_int_m + self.szy//2
            x_int[order-self.order_min, y_ix] = \
                np.interp(y_int_m, ypx[w_ix], xpx[w_ix])
            wave_int[order-self.order_min, y_ix] = \
                np.interp(y_int_m, ypx[w_ix], wave[w_ix])
            blaze_int[order-self.order_min, y_ix] = np.interp(
                y_int_m, ypx[w_ix], np.sinc(order_frac[w_ix])**2)
            # plt.plot(x_int[m-self.order_min,ix],y_int_m)
        # plt.axis( (-self.szx/2,self.szx/2,-self.szx/2,self.szx/2) )
        # plt.draw()
        ccd_centre = {'ccdx': ccdx, 'ccdy': ccdy, 'xpix_offset': xpix_offset,
                      'mean_v': mean_v}
        return x_int, wave_int, blaze_int, ccd_centre

    def spectral_format_with_matrix(self):
        """Create a spectral format, including a detector to slit matrix at
        every point.

        Returns
        -------
        x: (n_orders, n_y) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel
            and each order (m).
        w: (n_orders, n_y) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (n_orders, n_y) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        matrices: (n_orders, n_y, 2, 2) float array
            2x2 slit rotation matrices.
        """
        x_c, w_c, b_c, ccd_centre = self.spectral_format()
        x_xp, w_xp, dummy_0, dummy_1 = \
            self.spectral_format(xoff=-1e-3, ccd_centre=ccd_centre)
        x_yp, w_yp, dummy_0, dummy_1 = \
            self.spectral_format(yoff=-1e-3, ccd_centre=ccd_centre)
        dy_dyoff = np.zeros(x_c.shape)
        dy_dxoff = np.zeros(x_c.shape)
        # For the y coordinate, spectral_format output the wavelength at
        # fixed pixel, not the pixel at fixed wavelength. This means we need
        # to interpolate to find the slit to detector transform.
        isbad = w_c*w_xp*w_yp == 0
        for i in range(x_c.shape[0]):
            # w_ix = np.where(isbad[i, :] == False)[0]
            w_ix = np.where(np.logical_not(isbad[i, :]))[0]
            dy_dyoff[i, w_ix] = \
                np.interp(w_yp[i, w_ix], w_c[i, w_ix],
                          np.arange(len(w_ix))) - np.arange(len(w_ix))
            dy_dxoff[i, w_ix] = \
                np.interp(w_xp[i, w_ix], w_c[i, w_ix],
                          np.arange(len(w_ix))) - np.arange(len(w_ix))
            # Interpolation won't work beyond the end, so extrapolate manually
            # (why isn't this a numpy option???)
            dy_dyoff[i, w_ix[-1]] = dy_dyoff[i, w_ix[-2]]
            dy_dxoff[i, w_ix[-1]] = dy_dxoff[i, w_ix[-2]]

        # For dx, no interpolation is needed so the numerical derivative is
        # trivial...
        dx_dxoff = x_xp - x_c
        dx_dyoff = x_yp - x_c

        # flag bad data...
        x_c[isbad] = np.nan
        w_c[isbad] = np.nan
        b_c[isbad] = np.nan
        dy_dyoff[isbad] = np.nan
        dy_dxoff[isbad] = np.nan
        dx_dyoff[isbad] = np.nan
        dx_dxoff[isbad] = np.nan
        matrices = np.zeros((x_c.shape[0], x_c.shape[1], 2, 2))
        amat = np.zeros((2, 2))

        for i in range(x_c.shape[0]):
            for j in range(x_c.shape[1]):
                # Create a matrix where we map input angles to
                # output coordinates.
                amat[0, 0] = dx_dxoff[i, j]
                amat[0, 1] = dx_dyoff[i, j]
                amat[1, 0] = dy_dxoff[i, j]
                amat[1, 1] = dy_dyoff[i, j]
                # Apply an additional rotation matrix.
                # If the simulation was complete, this wouldn't be required.
                r_rad = np.radians(self.extra_rot)
                dy_frac = (j - x_c.shape[1]/2.0)/(x_c.shape[1]/2.0)
                extra_rot_mat = np.array([[np.cos(r_rad*dy_frac),
                                           np.sin(r_rad*dy_frac)],
                                          [-np.sin(r_rad*dy_frac),
                                           np.cos(r_rad*dy_frac)]])
                # print 'amat',amat
                # print 'extra_rot_mat',extra_rot_mat
                amat = np.dot(extra_rot_mat, amat)
                # We actually want the inverse of this
                # (mapping output coordinates back onto the slit).
                # print amat.shape
                if np.any(np.isnan(amat)):
                    #!!! Below this was zeros... !!!
                    matrices[i, j, :, :] = np.eye(len(amat))
                else:
                    matrices[i, j, :, :] = np.linalg.inv(amat)

        # Store these key spectral format parameters as object properties
        self.x = x_c
        self.wave = w_c
        self.blaze = b_c
        self.matrices = matrices
        self.stale_spectral_format = False
        return x_c, w_c, b_c, matrices

    def make_lenslets(self, fluxes=[], mode='', seeing=0.8, llet_offset=0, ifu=1):
        """Make an image of the lenslets with sub-pixel sampling.

        Parameters
        ----------
        fluxes: float array (optional)
            Flux in each lenslet

        mode: string (optional)
            'high' or 'std', i.e. the resolving power mode of the
            spectrograph.  Either mode or fluxes must be set.

        seeing: float (optional)
            If fluxes is not given, then the flux in each lenslet is
            defined by the seeing.

        llet_offset: int
            Offset in lenslets to apply to the input spectrum"""
        print("Computing a simulated slit image...")
        szx = self.im_slit_sz
        szy = 256
        fillfact = 0.98
        hex_scale = 1.15 #!!! This exists elsewhere !!!
        # equivalent to a 1 degree FWHM for an f/3 input ???
        # !!! Double-check !!!
        conv_fwhm = 30.0
        if len(fluxes) == N_HR_TOT:
            mode = 'high'
        elif len(fluxes) == N_SR_TOT:
            mode = 'std'
        elif len(mode) == 0:
            raise ValueError("Error: "+N_SR_TOT+" or "+N_HR_TOT+" lenslets needed... "
                             "or mode should be set")
        if mode == 'std':
            n_lenslets = N_SR_TOT
            lenslet_width = self.lenslet_std_size
            fibers = SRFibers(lenslet_width, self.microns_pix)
            xoffset, yoffset = fibers.xyoffsets(ifu)
        elif mode == 'high':
            n_lenslets = N_HR_TOT
            lenslet_width = self.lenslet_high_size
            fibers = HRFibers(lenslet_width, self.microns_pix)
            xoffset, yoffset = fibers.xyoffsets(ifu)
        else:
            raise ValueError("Error: mode must be std or high")

        # Some preliminaries...
        cutout_hw = int(lenslet_width/self.microns_pix*1.5)
        im_slit = np.zeros((szy, szx))
        x = np.arange(szx) - szx/2.0
        y = np.arange(szy) - szy/2.0
        xy_mesh = np.meshgrid(x, y)
        # radius and w_r enable the radius from the lenslet center to be indexed
        radius = np.sqrt(xy_mesh[0]**2 + xy_mesh[1]**2)
        w_r = np.where(radius < 2*lenslet_width/self.microns_pix)
        # g_frd is a Gaussian used for FRD
        g_frd = np.exp(-radius**2/2.0/(conv_fwhm/self.microns_pix/2.35)**2)
        g_frd = np.fft.fftshift(g_frd)
        g_frd /= np.sum(g_frd)
        gft = np.conj(np.fft.rfft2(g_frd))
        pix_size_slit = self.px_sz * (self.f_col / self.assym) \
            / self.f_cam * 1000.0 / self.microns_pix
        pix = np.zeros((szy, szx))
        pix[np.where((np.abs(xy_mesh[0]) < pix_size_slit/2) * \
                     (np.abs(xy_mesh[1]) < pix_size_slit/2))] = 1
        pix = np.fft.fftshift(pix)
        pix /= np.sum(pix)
        pix_ft = np.conj(np.fft.rfft2(pix))
        # Create some hexagons. We go via a "cutout" for efficiency.
        h_cutout = \
            optics.hexagon(
                szy, lenslet_width/self.microns_pix*fillfact/hex_scale)
        hbig_cutout = \
            optics.hexagon(szy, lenslet_width/self.microns_pix*fillfact)
        h_long = np.zeros((szy, szx))
        h_long_big = np.zeros((szy, szx))
        h_long[:, szx//2-szy//2:szx//2+szy//2] = h_cutout
        h_long_big[:, szx//2-szy//2:szx//2+szy//2] = hbig_cutout
        if len(fluxes) != 0:
            # If we're not simulating seeing, the image-plane is uniform,
            # and we only use the values of "fluxes" to scale the lenslet
            # fluxes.
            im_object = np.ones((szy, szx))
            # Set the offsets to zero because we may be simulating
            # e.g. a single Th/Ar lenslet and not starlight
            # (from the default xoffset etc)
            xoffset = np.zeros(len(fluxes), dtype=int)
            yoffset = np.zeros(len(fluxes), dtype=int)
        else:
            # If we're simulating seeing, create a Moffat function as our
            # input profile, but just make the lenslet fluxes uniform.
            im_object = np.zeros((szy, szx))
            im_cutout = optics.moffat2d(szy,  # noqa
                seeing * self.microns_arcsec / self.microns_pix / 2, beta=4.0)
            im_object[:, szx//2-szy//2:szx//2+szy//2] = im_cutout
            # Scale the image so the mean is 1.0, simply for rough numerical
            # consistency with the no-seeing case above.
            im_object /= im_object.mean()
            fluxes = np.ones(len(xoffset))

        # Go through the flux vector and fill in each lenslet.
        cutoutx = [max(0, szx//2 - cutout_hw), min(szx//2 + cutout_hw, szx)]
        cutouty = [max(0, szy//2 - cutout_hw), min(szy//2 + cutout_hw, szy)]
        for i, flux in enumerate(fluxes):
            im_one = np.zeros((szy, szx))
            im_cutout = np.roll(np.roll(im_object, yoffset[i], axis=0),
                                xoffset[i], axis=1) * h_long
            im_cutout = im_cutout[cutouty[0]:cutouty[1],
                                  cutoutx[0]:cutoutx[1]]
            prof = optics.azimuthal_average(im_cutout, returnradii=True,
                                            binsize=1)
            prof = (prof[0], prof[1] * flux)
            xprof = np.append(np.append(0, prof[0]), np.max(prof[0])*2)
            yprof = np.append(np.append(prof[1][0], prof[1]), 0)
            im_one[w_r] = np.interp(radius[w_r], xprof, yprof)
            im_one = np.fft.irfft2(np.fft.rfft2(im_one)*gft)*h_long_big
            im_one = np.fft.irfft2(np.fft.rfft2(im_one)*pix_ft)
            # !!! The line below could add tilt offsets...
            # important for PRV simulation !!!
            # im_one = np.roll(np.roll(im_one, tilt_offsets[0,i], axis=1),\
            #tilt_offsets[1,i], axis=0)*h_long_big
            the_shift = \
                int((llet_offset + i - n_lenslets/2.0) *
                    lenslet_width/self.microns_pix)
            im_slit += np.roll(im_one, the_shift, axis=1)

        # Now normalise the slit image to have a sum of 1.0. This is the easiest
        # type of normalisation to deal with for conversion to physical units.
        im_slit /= np.sum(im_slit)
        return im_slit

    def get_solar_spectrum(self, flux_scale=100.0):
        """Extract a default solar spectrum from a file and return as a 2D array"""
        data = pf.getdata(os.path.join(LOCAL_DIR, 'data/ardata.fits.gz'))
        spectrum = np.array([np.append(0.35, data['WAVELENGTH']) / 1e4,
                             np.append(0.1, np.maximum(data['SOLARFLUX'],0))*flux_scale])
        return spectrum


    def simulate_image(self, x, wave, blaze, matrices, im_slit, spectrum=None,
                       n_x=0, xshift=0.0, yshift=0.0, radvel=0.0, return_check=False):
        """Simulate a spectrum on the CCD.

        Parameters
        ----------
        x, wave, blaze, matrices: float arrays
            See the output of spectral_format_with_matrix
        im_slit: float array
            See the output of make_lenslets
        spectrum: (2,nwave) array (optional)
            An input spectrum, arbitrarily gridded (but at a finer
            resolution than the spectrograph resolving power).
            If not given, a solar spectrum is used.
        n_x: float
            Number of x (along-slit) direction pixels in the image.
            If not given or zero, a square CCD is assumed.
        xshift: float
            Bulk shift to put in to the spectrum along the slit.
        yshift: floatDROP SCE
            NOT IMPLEMENTED
        radvel: float
            Radial velocity in m/s.
        return_check: bool
            Also return the input spectra, interpolated onto the pixel
            grid for each order.
        """
        # If no input spectrum, use the sun.
        if (spectrum is None) or len(spectrum) == 0:
            spectrum = self.get_solar_spectrum()

        n_orders = x.shape[0]
        n_y = x.shape[1]
        if n_x == 0:
            n_x = n_y
        image = np.zeros((n_y, n_x))
        if return_check:
            check_image = np.zeros((n_orders, n_y))
        # Simulate the slit image within a small cutout region.
        cutout_xy = np.meshgrid(np.arange(81)-40, np.arange(7)-3)
        # Loop over orders. Use a status line where the numbers change but new
        # lines aren't created.
        for i in range(n_orders):
            for j in range(n_y):
                if x[i, j] != x[i, j]:
                    continue
                # We are looping through y pixel and order.
                # The x-pixel is therefore non-integer.
                # Allow an arbitrary shift of this image.
                the_x = x[i, j] + xshift
                # Create an (x,y) index of the actual pixels we want
                # to index.
                cutout_shifted = (cutout_xy[0].copy() + int(the_x) + n_x//2,
                                  cutout_xy[1].copy() + j)
                w_ix = np.where((cutout_shifted[0] >= 0) *
                                (cutout_shifted[1] >= 0) *
                                (cutout_shifted[0] < n_x) *
                                (cutout_shifted[1] < n_y))
                cutout_shifted = (cutout_shifted[0][w_ix],
                                  cutout_shifted[1][w_ix])
                flux = np.interp(wave[i, j]*(1 + radvel/299792458.0),
                                 spectrum[0], spectrum[1],
                                 left=0, right=0)
                if return_check:
                    check_image[i, j] = flux
                # Rounded to the nearest microns_pix, find the co-ordinate
                # in the simulated slit image corresponding to each pixel.
                # The co-ordinate order in the matrix is (x,y).
                xy_scaled = np.dot(
                    matrices[i, j],
                    np.array([cutout_xy[0][w_ix] + int(the_x) - the_x,
                              cutout_xy[1][w_ix]]) / self.microns_pix
                ).astype(int)
                slit_y = xy_scaled[1] + im_slit.shape[0]//2
                slit_x = xy_scaled[0] + im_slit.shape[1]//2
                w_ix = np.where((slit_x >= 0) * (slit_y >= 0) *
                                (slit_x < im_slit.shape[1]) *
                                (slit_y < im_slit.shape[0]))
                cutout_shifted = (cutout_shifted[0][w_ix],
                                  cutout_shifted[1][w_ix])
                slit_x = slit_x[w_ix]
                slit_y = slit_y[w_ix]
                # Find the flux scaling. The slit image is normalised to 1, and we
                # normalize the re-sampled slit image to 1 by multiplying by the
                # determinant of the slit-to-detector matrix.
                image[cutout_shifted[1], cutout_shifted[0]] += blaze[i, j] * \
                    flux * im_slit[slit_y, slit_x] * np.linalg.det(matrices[i,j])
            print('\r Doing orders {0:d} to {1:d}. Done order {2:d}'.format(self.order_min, self.order_max, i+self.order_min), end='\r')
            sys.stdout.flush()
            #time.sleep(0.1)
            #outline.write('Doing orders {0:d} to {1:d}. Done order {2:d}\r'.format(self.order_min, self.order_max, i+self.order_min))
            #outline.flush()
        print('\n')
        #outline.write('\n')
        #Old code was: print('Done order: {0}'.format(i + self.order_min))
        if return_check:
            return (image, check_image)
        else:
            return image

    def sky_background(self, mode):
        """Calculate a sky spectrum.
        Parameters
        ----------
        mode: string
            Either 'std' or 'high' depending on the IFU mode in use.

        Returns
        -------
        wave, flux: sky spectrum (wavelength in um, flux in photons/s)
        """

        # Input checks
        if mode not in self.RES_OPTIONS:
            raise ValueError('Mode must be one of %s' % (
                ', '.join(self.RES_OPTIONS),
            ))

        # Load the UVES sky
        # Flux is in 1e-16 erg / (s*A*cm^2*arcs^2)
        bgwave, bgflux = load_sky_from_dir(os.path.join(LOCAL_DIR, 'data'))

        # Convert flux to phot/s/A/cm^2/arcsec^2
        bgflux *= 1e-16
        bgflux /= PLANCK_H*LIGHT_C/bgwave

        # Calculate area per fiber in arcsec^2
        # This is the area of a hexagon on the focal plane, where the scale
        # is 610 um/arcsec.
        if mode == 'high':
            fiber_area = S32 * 240.0**2 / (610.0 ** 2)
        elif mode == 'std':
            fiber_area = S32 * 144.0**2 / (610.0 ** 2)

        # Convert phot/s/A/cm^2/arcsec^2 into phot/s/A/cm^2
        bgflux *= fiber_area

        # Convert phot/s/A/cm^2 to phot/s/A
        # (assuming an 8m mirror)
        bgflux *= math.pi * 400**2

        # Convert to phot/s
        bgflux = bgflux[:-1] * np.diff(bgwave)

        # Convert wavelength to microns
        bgwave = bgwave[:-1] / 1e4

        # plt.plot(bgwave, bgflux)
        # plt.show()

        # Return flux and wavelength
        return np.array([bgwave, bgflux])

    def simulate_frequency_noise(self, freq, mean, std, binmode=(1, 1),
                                 overscan=0):
        """ Simulate an image with noise of specific frequencies in it. """
        shape = (self.szx//binmode[1], (self.szy+2*overscan)//binmode[0])
        return frequency_noise(
            freq, self.sample_rate, shape, mean=mean, std=std)

    def simulate_gradient(self, theta, mean, std, binmode=(1, 1), overscan=0):
        """ Simulate an image with a gradient across it, at angle theta.
            The image has mean=0 and std=1. """
        x, y = np.meshgrid(np.arange((
            self.szy+2*overscan)/binmode[0]), np.arange(self.szx/binmode[1]))
        image = x * np.sin(theta) + y * np.cos(theta)
        image = (image - image.mean()) / image.std()
        return image * std + mean

    def simulate_flatfield(self, mean, std, overscan=0):
        """ Simulate a flatfield. """
        shape = (self.szx, self.szy+2*overscan)
        flat = np.random.normal(loc=mean, scale=std, size=shape)
        idx = tuple([np.random.choice(dm, 300) for dm in shape])
        flat[idx] = 0.0  # simulate dead pixels
        return flat

    def blank_frame(self, namps=[1, 1]):
        """ Return an entirely blank frame, of the correct size """
        image = np.zeros((self.szx, self.szy))
        images = split_image(image, namps)
        return images

    def simulate_frame(self, duration=0.0, output_prefix='test_',
                       spectrum_in=None, xshift=0.0, yshift=0.0, radvel=0.0,
                       rv_thar=0.0, flux=1, rnoise=3.0, gain=[1.0],
                       bias_level=10, overscan=0, namps=[1, 1], use_thar=True,
                       res='high', add_sky=True, return_image=False,
                       thar_flatlamp=False, flatlamp=False, obstype=None,
                       objname=None, additive_noise={}, scaling=None,
                       seeing=0.8, data_label=1,
                       utstart=datetime.datetime.utcnow(), binning=(1, 1)):
        """Simulate a single frame.

        TODO (these can be implemented manually using the other functions):
        1) Variable seeing (the slit profile is currently fixed)
        2) Standard resolution mode.
        3) Sky
        4) Arbitrary input spectra

        Parameters
        ----------
        duration: float (optional)
            Duration of the exposure in seconds

        output_prefix: string (optional)
            Prefix for the output filename.

        spectrum_in: array of shape (2,n) where n is the number of
            spectral samples input spectrum.
            spectrum[0] is the wavelength array, spectrum[1] is the flux array,
            in units of recorded photons/spectral pixel/s, taking into account
            the full instrument at order centers (i.e. excluding the blaze
            function)

        xshift: float (optional)
            x-direction (along-slit) shift.

        yshift: float (optional)
            y-direction (spectral direction) shift.

        radvel: float (optional)
            Radial velocity in m/s for the target star with respect to
            the observer.

        rv_thar: float (optional)
            Radial velocity in m/s applied to the Thorium/Argon source.  It is
            unlikely that this is useful (use yshift instead for common shifts
            in the dispersion direction).

        flux: float (optional) Flux multiplier for the reference spectrum to
            give photons/pix/s.

        rnoise: float (optional)
            Readout noise in electrons/pix

        gain: float[namps] (optional)
            Gain in electrons per ADU, per amplifier (or a single scalar for all
            amps).

        bias_level: float (optional)
            Bias level in electrons, per amplifier (or a single scalar for
            all amps).

        overscan: int (optional)
            number of columns of overscan

        namps: int[2] (optional)
            number of readout amps in each direction

        use_thar: bool (optional)
            Is the Thorium/Argon lamp in use?

        obstype: string
            The value for the OBSTYPE FITS keyword

        obsname: string
            The name of the target to be created. Will populate the OBJECT
            keyword.
        
        res: string (optional)
            Can be 'high' or 'std' for the resolution mode.

        add_sky: bool (optional)
            Is the sky background added to the frame?

        return_image: bool (optional)
            Do we return an image as an array? The fits file is always written.

        data_label: int (optional)
            Indicates which frame of a sequence this is; written to the DATALAB
            keyword

        utstart: datetime (optional)
            The UT start time from which to generate all time-related header
            keywords

        binning: 2-tuple of ints (optional)
            The binning mode to apply to the output data frame
        """

        # Input checks
        if res not in self.RES_OPTIONS:
            raise ValueError('Mode must be one of %s' % (
                ', '.join(self.RES_OPTIONS),
            ))

        # If no input spectrum, use the sun with default scaling.
        if (spectrum_in is None) or len(spectrum_in) == 0:
            spectrum_in = self.get_solar_spectrum()

        # Scale the spectrum as required.
        spectrum = spectrum_in * flux

        # Compute the spectral format if we need to
        if self.stale_spectral_format:
            x, wave, blaze, matrices = self.spectral_format_with_matrix()
        else:
            x = self.x
            wave = self.wave
            blaze = self.blaze
            matrices = self.matrices

        # Deal with the values that can be scalars or arrays.
        # Turn them into arrays of the right size if they're scalars.
        if np.isscalar(gain):
            gain = np.ones(namps[0]*namps[1]) * gain

        if np.isscalar(bias_level):
            bias_level = np.ones(namps[0]*namps[1]) * bias_level

        if np.isscalar(rnoise):
            rnoise = np.ones(namps[0]*namps[1]) * rnoise

        if scaling is None:
            scaling = np.ones(namps[0]*namps[1])
        elif np.isscalar(scaling):
            scaling = np.ones(namps[0]*namps[1]) * scaling

        # Make sure duration is a float so that it doesn't mess
        # up our calculations later on
        duration = float(duration)

        hdr = pf.Header()
        hdr['TARGET1'] = 0  # ifu1 points to nothing in particular to begin with
        hdr['TARGET2'] = 0  # ifu2 points to nothing in particular to begin with

        # WARNING: Logic below can be neatened... (e.g. doubled-up lines of
        # code)
        if res == 'high':
            slit_fluxes = np.ones(N_HR_SCI+N_HR_SKY) if flatlamp else []
            im_slit = self.make_lenslets(
                fluxes=slit_fluxes, mode=res, seeing=seeing, llet_offset=2)
            if obstype == 'OBJECT' or obstype == 'STANDARD':
                hdr['TARGET1'] = 2  # ifu1 points to an object
            image, check_image = self.simulate_image(
                x, wave, blaze, matrices, im_slit, spectrum=spectrum,
                n_x=self.szx, xshift=xshift, radvel=radvel, return_check=True)

            # Create the slit viewing camera image for the science flux.
            if self.slitv is not None:
                self.slitv.create_slitview_frames(
                    im_slit, spectrum, self.slitview_wave, self.slitview_frac,
                    self.slitview_offset, res)

            if use_thar:
                # Create an appropriately convolved Thorium-Argon spectrum
                thar_spect = thar_spectrum()
                '''
                plt.clf()
                plt.plot(thar_spect[0], thar_spect[1])
                plt.show()
                '''
                if thar_flatlamp:
                    thar_spect[1][:] = 10

                # Now that we have our spectrum, create the Th/Ar image.
                slit_fluxes = np.ones(1)
                im_slit2 = self.make_lenslets(
                    fluxes=slit_fluxes, mode=res, llet_offset=0)
                image += self.simulate_image(
                    x, wave, blaze, matrices, im_slit2, spectrum=thar_spect,
                    n_x=self.szx, xshift=xshift, radvel=rv_thar)

                # Add in the Th/Ar spectrum
                if self.slitv is not None:
                    self.slitv.create_slitview_frames(
                        im_slit2, thar_spect, self.slitview_wave,
                        self.slitview_frac, self.slitview_offset, res)

        elif res == 'std':
            if flatlamp:
                slit_fluxes = np.ones(N_SR_TOT)
                im_slit = self.make_lenslets(
                    fluxes=slit_fluxes, mode=res, seeing=seeing, llet_offset=0)

            else:
                im_slit = self.make_lenslets(
                    fluxes=[], mode=res, seeing=seeing, llet_offset=0, ifu=1)
                if obstype == 'OBJECT' or obstype == 'STANDARD':
                    hdr['TARGET1'] = 2  # ifu1 points to an object
                if self.dual_target:
                    im_slit += self.make_lenslets(
                        fluxes=[], mode=res, seeing=seeing, llet_offset=N_SR_SCI +
                        N_SR_SKY, ifu=2)
                    if obstype == 'OBJECT' or obstype == 'STANDARD':
                        hdr['TARGET2'] = 2  # ifu2 points to an object: shouldn't
                                            # really be the same one, but oh well

            image, check_image = self.simulate_image(
                x, wave, blaze, matrices, im_slit, spectrum=spectrum,
                n_x=self.szx, xshift=xshift, radvel=radvel, return_check=True)

            # Create the slit viewing camera image for the science flux.
            if self.slitv is not None:
                self.slitv.create_slitview_frames(
                    im_slit, spectrum, self.slitview_wave, self.slitview_frac,
                    self.slitview_offset, res)

        else:
            # mode input check above makes this test unecessary/redundant
            print("ERROR: unknown resolution.")
            raise UserWarning

        # Prevent any interpolation errors (negative flux)
        # prior to adding noise.
        image = np.maximum(image, 0)

        # Scale to photons
        image = np.random.poisson(duration * image).astype(float)

        if add_sky and obstype != 'BIAS' and obstype != 'DARK':
            # Calculate the sky spectrum - the flux we calculate is per fiber.
            sky_spect = self.sky_background(res)
            # Put it into all the fibers equally
            if res == 'high':
                slit_fluxes = np.ones(N_HR_SCI + N_HR_SKY)
                im_slit2 = self.make_lenslets(
                    fluxes=slit_fluxes, mode=res, llet_offset=2)
                hdr['TARGET1'] = 1 if hdr['TARGET1'] == 0 else hdr['TARGET1']

            elif res == 'std':
                slit_fluxes = np.ones(N_SR_TOT)
                im_slit2 = self.make_lenslets(
                    fluxes=slit_fluxes, mode=res, llet_offset=0)
                hdr['TARGET1'] = 1 if hdr['TARGET1'] == 0 else hdr['TARGET1']
                hdr['TARGET2'] = 1 if hdr['TARGET2'] == 0 else hdr['TARGET2']

            sky_image = self.simulate_image(
                x, wave, blaze, matrices, im_slit2, spectrum=sky_spect,
                n_x=self.szx, xshift=xshift, radvel=0.0)

            sky_image = np.maximum(0, sky_image)
            # Add to the science image
            image += np.random.poisson(duration * sky_image)
            # And to the slit viewing image
            if self.slitv is not None:
                self.slitv.create_slitview_frames(
                    im_slit2, sky_spect, self.slitview_wave, self.slitview_frac,
                    self.slitview_offset, res)

        # We ignore atmospheric transmission and instrument throughput in this
        # routine, because the input spectrum is assumed to include all of these
        # effects other than the blaze function.

        if self.cosmics:
            # Add cosmic rays
            # We hard code 2.0 rays/s/cm^2
            # A shield that only allows rays that originate within 10 degrees
            # of the zenith.

            # Pixel size of 15 um x 15 um x pixel_depth um
            cosmic_img = cosmic.cosmic(image.shape,     # Image shape
                                       duration,        # Exposure length
                                       10,              # CR shield angle
                                       2.0,             # CR/s/cm/cm
                                       False,           # Use effect area mask
                                       [15, 15, self.thick]  # Pxl sz (microns)
                                       )
            image += cosmic_img
            # import pdb; pdb.set_trace() #!!!MJI!!!
            # no_cr_pix = np.count_nonzero(cosmic_img)
        else:
            cosmic_img = np.zeros(image.shape)

        # Compute dark current (3 e/pix/hour)
        dcurrent = np.float64(np.random.poisson(
            np.ones_like(image) * 3.0 * duration/3600.0))

        if self.hpplane:
            hp_img = self.hot_pix.toarray().T
            hp_img = to_ushort(hp_img)
            hphdu = pf.HDUList(pf.PrimaryHDU(data=hp_img, header=pf.Header()))
            hpims = split_image(hp_img, namps)
            hpims = [add_overscan(i, overscan) for i in hpims]
            hpims = [to_ushort(i) for i in hpims]
            for hpim in hpims:
                hphdu.append(pf.ImageHDU(data=hpim, header=pf.Header()))
            hphdu.writeto(self.arm + '_HP.fits', overwrite=True)
            self.hpplane = False

        # Add hot pixels to the dark current image
        dcurrent += self.hot_pix.toarray() * np.std(dcurrent)

        # Add dark current and hot pixels to frame
        image += dcurrent
        # Or, have a ludicrously high dark current (1500 e/pix/hr)
        # image += np.random.poisson(
        #     np.ones_like(image) * 1500.0 * duration / 3600.0)
        # Or, have a high dark current on half the CCD and zero on the other
        # half
        # hlfw = image.shape[-1] / 2
        # image += np.concatenate(
        #     (
        #         np.random.poisson(np.ones_like(image[:, :hlfw]) * 36000.0 *
        #                           duration / 3600.0),
        #         np.zeros_like(image[:, hlfw:]),
        #     ), axis=-1)

        # MCW 190814 - Flip the image y-axis (axis 0)
        image = image.T[::-1, :]  # transpose image for conventional axes
        # Split the image into sections for each readout amplifier
        ampims, cxl, cxh, cyl, cyh = split_image(
            image, namps, return_headers=True)

        # MCW 190814 - Flip the image y-axis (axis 0)
        cosmic_img = cosmic_img.T[::-1, :]
        cosims = split_image(cosmic_img, namps)
        cosims = [add_overscan(i, overscan) for i in cosims]
        if self.split and self.crplane and self.cosmics:
            cosims = [to_ushort(i) for i in cosims]
            crhdu = pf.HDUList(pf.PrimaryHDU(header=pf.Header()))
            for i, cosim in enumerate(cosims):
                crhdr = pf.Header()
                crhdr['DETSIZE'] = ("[1:%d,1:%d]" % (
                    cosim.shape[1], cosim.shape[0]), 'Detector size')
                detsec = "[%d:%d,%d:%d]" % (cxl[i]+1, cxh[i], cyl[i]+1, cyh[i])
                crhdr['DETSEC'] = (detsec, 'Detector section(s)')
                cosdsec = "[1:%d,1:%d]" % (
                    cosim.shape[1]-overscan, cosim.shape[0])
                crhdr['DATASEC'] = (cosdsec, 'Data section(s)')
                crhdr['TRIMSEC'] = (cosdsec, 'Trim section(s)')
                if overscan > 0: crhdr['BIASSEC'] = ("[%d:%d,1:%d]" % (
                    cosim.shape[1]-overscan+1, cosim.shape[1],
                    cosim.shape[0]), 'Bias section(s)')
                crhdu.append(pf.ImageHDU(data=cosim, header=crhdr, name='SCI'))

            print('Writing ' + output_prefix + self.arm + '_CR.fits')
            crhdu.writeto(output_prefix + self.arm + '_CR.fits', overwrite=True)

        # Write the check frame if required
        if self.split and self.check and obstype != 'BIAS' and obstype != 'DARK':
            chkhdu = pf.HDUList(pf.PrimaryHDU(data=check_image))
            chkhdu.writeto(output_prefix + self.arm + '_chk.fits', overwrite=True)

        # some datetimes for use in creating the header below
        ltnow = utstart.replace(tzinfo=tz.tzutc())
        ltnow = ltnow.astimezone(tz.tzlocal())

        # Now create our fits image!
        # By adding the DETSIZE and DETSEC keywords we can open the
        # raw image in ds9 as an IRAF mosaic.
        hdr['OBSERVAT'] = (
            'Gemini-South', 'Name of telescope (Gemini-North|Gemini-South)')
        hdr['TELESCOP'] = 'Gemini-South'
        hdr['INSTRUME'] = ('GHOST', 'Instrument used to acquire data')
        if obstype=='STANDARD':
            hdr['OBSTYPE'] = ('OBJECT', 'Observation type')
        else:
            hdr['OBSTYPE'] = (obstype, 'Observation type')

        # required by the calibration manager
        hdr['RAWPIREQ'] = 'yes'  # 'no', 'check', 'unknown'
        hdr['RAWGEMQA'] = 'usable'  # 'bad', 'check', 'unknown'

        # these are removed from the primary header when writing the MEF;
        # conversely the MEF splitter primitive will (may?) need to add
        # these to the split individual frame files in order for them to be
        # correctly categorized by the type system
        detsz = "[1:%d,1:%d]" % (image.shape[1], image.shape[0])
        hdr['DETSIZE'] = (detsz, 'Detector size')
        hdr['DETTYPE'] = (self.dettype, 'Detector array type')
        hdr['CAMERA'] = (self.arm.upper(), 'Camera name')

        if objname==None:
            # populate OBJECT keyword
            obj = dict(FLAT='GCALflat', ARC='ThAr', OBJECT='V492 Car',
                       BIAS='Bias', DARK='Dark', SKY='')  # noqa
            hdr['OBJECT'] = (obj[obstype], 'Object Name')
        else:
            hdr['OBJECT'] = (objname, 'Object Name')
        # populate OBSCLASS keyword
        obsclass = dict(FLAT='partnerCal', ARC='partnerCal',  # noqa
                        OBJECT='science', BIAS='dayCal', DARK='dayCal',
                        SKY='', STANDARD='partnerCal')
        hdr['OBSCLASS'] = (obsclass[obstype], 'Observe class')

        # populate GEMPRGID, OBSID, and DATALAB keywords
        prgid = 'GS-2016B-Q-20'  # just choose a base at random
        hdr['GEMPRGID'] = (prgid, 'Gemini programme ID')
        obsid = dict(
            high=dict(BIAS='1', DARK='5', FLAT='4', ARC='10', OBJECT='9',
                      SKY='7', STANDARD='9'),  # noqa
            std=dict(BIAS='1', DARK='5', FLAT='3', ARC='2', OBJECT='8',
                     SKY='6', STANDARD='8'),  # noqa
            )
        hdr['OBSID'] = (
            prgid+'-'+obsid[res][obstype], 'Observation ID / Data label')
        hdr['DATALAB'] = (prgid+'-'+obsid[res][obstype]+'-' + (
            '%03d' % data_label), 'DHS data label')

        # keywords that vary by time and position: first start with
        # a spot in the sky that's always visible from Gemini South
        my_ra, my_dec = 149.2442500, -69.1009167  # V492 Car
        spot = SkyCoord(my_ra, my_dec, unit=u.deg, frame='fk5')
        # location and UTC time at Gemini South
        gs_locn = EarthLocation(
            lat=-30.24075*u.deg, lon=-70.736693*u.deg, height=2722*u.m)
        gs_time = Time(ltnow, location=gs_locn) - 3*u.hour
        aa_start = AltAz(obstime=gs_time, location=gs_locn)
        aa_end = AltAz(obstime=gs_time+duration*u.second, location=gs_locn)
        horz_start = spot.transform_to(aa_start)
        horz_end = spot.transform_to(aa_end)
        hdr['RA'] = (my_ra, 'Right Ascension')
        hdr['DEC'] = (my_dec, 'Declination of Target')
        hdr['ELEVATIO'] = (horz_start.alt.deg, 'Current Elevation')
        hdr['AZIMUTH'] = (horz_start.az.deg, 'Current Azimuth')
        # TODO: use apparent sidereal time below instead of mean?
        hra = gs_time.sidereal_time('mean') - my_ra*u.deg  # hour angle
        hdr['HA'] = (hra.to_string(unit=u.deg, sep=':'), 'Telescope hour angle')
        hdr['AMSTART'] = (horz_start.secz.value, 'Airmass at start of exposure')
        hdr['AMEND'] = (horz_end.secz.value, 'Airmass at end of exposure')
        # assumes linear transition from AMSTART to AMEND which is likely wrong
        mean_airmass = (horz_end.secz.value + horz_start.secz.value)/2
        hdr['AIRMASS'] = (mean_airmass, 'Mean airmass for the observation')

        # keywords that are constant
        hdr['EQUINOX'] = (2000., 'Equinox of coordinate system')
        hdr['HUMIDITY'] = (39., 'The Relative Humidity (fraction, 0..101)')
        hdr['TAMBIENT'] = (8.8, 'The ambient temp (C)')
        hdr['PRESSURE'] = (546.95412, 'The atmospheric pressure (mm Hg)')
        hdr['PA'] = (90., 'Sky Position Angle at start of exposure')
        hdr['IAA'] = (359.78, 'Instrument Alignment Angle')
        hdr['CRPA'] = (-187.344883172036, 'Current Cass Rotator Position Angle')
        # hdr['CRFOLLOW'] = (0, '')  # TODO: where can I get this info?

        # time-related keywords
        hdr['DATE-OBS'] = (
            utstart.strftime("%Y-%m-%d"), 'UT date at observation start')
        hdr['UTSTART'] = (utstart.strftime("%H:%M:%S.%f")[:-3],
            'UT time at observation start')  # noqa
        hdr['UTEND'] = (
            (utstart + datetime.timedelta(seconds=duration))
            .strftime("%H:%M:%S.%f")[:-3], 'UT time at observation end')
        hdr['LT'] = (ltnow.strftime("%H:%M:%S.%f")[:-3],
            'Local time at start of observation')  # noqa

        # resolution-related keywords
        if obstype != 'BIAS' and obstype != 'DARK':
            hdr['SMPNAME'] = ('HI_ONLY' if res == 'high' else 'LO_ONLY')
            hdr['SMPPOS'] = (1 if res == 'high' else 2)

        hdulist = pf.HDUList(pf.PrimaryHDU(header=hdr))

        binmodes = [binning]
        if self.split and obstype not in ['DARK', 'FLAT', 'ARC']:
            binmodes = [(1, 1), (1, 2), (1, 8), (2, 4), (2, 8)]
        elif obstype == 'BIAS' and binning != (1, 1):
            binmodes.insert(0, (1, 1))

        for expid, binmode in enumerate(binmodes, start=1):
            opims = deepcopy(ampims)

            # slap on an overscan region
            opims = [add_overscan(i, overscan) for i in opims]

            # scale each amp by its scaling factor (the 0-valued overscan region
            # will remain 0)
            opims = [i*s for i, s in zip(opims, scaling)]

            # apply binning
            opims = [apply_binning(i, binmode) for i in opims]
            oscan = overscan / binmode[0]  # the new binned overscan width

            # Add read noise (in electrons) for each amplifier
            opims = [i+r*np.random.normal(size=i.shape)
                for i, r in zip(opims, rnoise)]  # noqa

            # convert electrons to ADU by dividing by the gain (in e/ADU)
            opims = [i/g for i, g in zip(opims, gain)]

            # Add in the additive noise (this is assumed to be electronic noise,
            # so it's in ADUs)
            if additive_noise and binmode in additive_noise:
                opims = [i+n for i, n in zip(opims, additive_noise[binmode])]

            # add in each amp's bias level
            opims = [i+b for i, b in zip(opims, bias_level)]

            # FIXME Apply non-linearity

            # Convert to unsigned short, and deal with saturation
            opims = [to_ushort(i) for i in opims]

            if return_image:
                return opims

            for i, ampim in enumerate(opims):
                xhdr = pf.Header()
                xhdr['EXTNAME'] = ('SCI', 'extension name')
                xhdr['BUNIT'] = 'ADU'
                # xhdr['BSCALE'] = 1  # added by default
                # xhdr['BZERO'] = 32768  # added by default
                xhdr['CAMERA'] = (self.arm.upper(), 'Camera name')
                xhdr['RDNOISE'] = (rnoise[i], 'Readout noise')
                xhdr['GAIN'] = (gain[i], 'Amplifier gain')
                xhdr['EXPID'] = (expid, 'Exposure ID')
                xhdr['PCOUNT'] = (0, 'Required keyword; must = 0')
                xhdr['GCOUNT'] = (1, 'Required keyword; must = 1')

                xhdr['EXPUTST'] = (utstart.strftime("%H:%M:%S.%f")[:-3],
                    'UT time at exposure start')  # noqa
                xhdr['EXPUTEND'] = (
                    (utstart + datetime.timedelta(seconds=duration))
                    .strftime("%H:%M:%S.%f")[:-3], 'UT time at exposure end')
                # needed by gemcombine.cl (used by stackFrames())
                xhdr['EXPTIME'] = (duration, 'Exposure time in seconds')

                xhdr['DETSIZE'] = (detsz, 'Detector size')
                xhdr['CCDSIZE'] = (detsz, 'CCD size')

                xhdr['CCDNAME'] = (self.dettype, 'CCD name')
                xhdr['AMPNAME'] = ('ABCD'[i], 'Amplifier name')
                xhdr['AMPSIZE'] = ("[1:%d,1:%d]" % (
                    ampim.shape[1], ampim.shape[0]), 'Amplifier size')

                detsec = "[%d:%d,%d:%d]" % (cxl[i]+1, cxh[i], cyl[i]+1, cyh[i])
                xhdr['DETSEC'] = (detsec, 'Detector section(s)')
                xhdr['CCDSEC'] = (detsec, 'CCD section(s)')

                xhdr['CCDSUM'] = (
                    ' '.join(map(str, binmode)), 'CCD pixel summing')
                # TODO: what's a realistic value for readout time?
                xhdr['DARKTIME'] = (duration + 5., 'Dark time (seconds)')

                datasec = "[1:%d,1:%d]" % (ampim.shape[1]-oscan, ampim.shape[0])
                xhdr['DATASEC'] = (datasec, 'Data section(s)')
                xhdr['TRIMSEC'] = (datasec, 'Trim section(s)')

                # TODO: where can I find these values?
                xhdr['DSPTIMBN'] = ('???', 'ARC timing board dsp code')
                xhdr['DSPTIMBV'] = (0, 'ARC timing board dsp version')

                if oscan > 0:
                    xhdr['BIASSEC'] = ("[%d:%d,1:%d]" % (
                        ampim.shape[1]-oscan+1, ampim.shape[1],
                        ampim.shape[0]), 'Bias section(s)')

                if self.cosmics:
                    xhdr['NCRPIX'] = np.count_nonzero(cosims[i])

                hdulist.append(pf.ImageHDU(data=ampim, header=xhdr, name='SCI'))

            if self.split:
                bins = 'x'.join(map(str, binmode))
                newfilename = output_prefix + bins + '_' + self.arm + '.fits'
                print('Writing ' + newfilename)
                hdulist.writeto(newfilename, overwrite=True)
                # clear hdulist in preparation for next binning mode
                hdulist = pf.HDUList(pf.PrimaryHDU(header=hdr))
                continue

        return None if self.split else hdulist


class Ghost(object):
    """A class to encapsulate both Arm objects that represent the two arms of
    the spectrograph. The initialisation function takes the fixed parameters for
    the two detectors (the things that won't change between exposures).

    Parameters
    ----------
    rnoise: float
        The read noise for both detectors (assumed to be the same for now)

    gain: float
        The gain for both detectors (assumed to be the same for now)

    namps: int[2]
        The number of amplifiers in the row and column directions for both
        detectors (assumed to be the same)

    overscan: int
        The number of columns of overscan for each amplifier (assumed to be
        the same for all amplifiers on both detectors)

    split: bool
        Indicate whether to produce a single MEF containing the frames taken
        during the observation, or to generate individual files per frame

    bias_level: float
        The bias level for both detectors (assumed to be the same for now)

    additive_noise: dict{'red'/'blue': float[namps, ampsize]}
        Additive noise (in ADUs) that will be added to every exposure.

    scaling: dict{'red'/'blue': float[namps, ampsize]}
        The scaling factor (i.e. pixel-to-pixel variation) for every exposure

    cosmics: bool
        Whether to add cosmic rays to frames or not

    crplane: bool
        Output a fits file containing locations where CRs were injected?

    hpplane: bool
        Output a fits file containing locations of hot pixels?

    split: bool
        Indicate whether to produce a single MEF containing the frames taken
        during the observation, or to generate individual files per frame

    check: bool
        Output a fits file containing input spectra interpolated onto
        the pixel grid for each order?

    dual_target: bool
        Output STD res fits files that contain object spectra in both IFUs
        (Ignored for HIGH res mode)
    """

    def __init__(self, rnoise, gain, namps, overscan, bias_level,
                 additive_noise, scaling, cosmics=True, crplane=False,
                 hpplane=False, split=False, check=False, dual_target=False):

        self.rnoise = rnoise
        self.gain = gain
        self.namps = namps
        self.overscan = overscan
        self.bias_level = bias_level
        self.additive_noise = additive_noise
        self.scaling = scaling
        self.cosmics = cosmics
        self.crplane = crplane
        self.hpplane = hpplane
        self.split = split
        self.check = check

        self.slitv = SlitViewer(cosmics, crplane, split)
        common = dict(
            slit_viewer=self.slitv, cosmics=cosmics, crplane=crplane,
            hpplane=hpplane, split=split, check=check, dual_target=dual_target)
        # pylint: disable=star-args
        self.blue = Arm('blue', **common)
        self.red = Arm('red', **common)


    def get_standard_spectrum(self, std='hd160617.fits'):
        """Extract the spectrum of one of the standard stars in the data folder"""
        try:
            data = pf.getdata(os.path.join(LOCAL_DIR, 'data/standards/'+std))
        except:
            return 'Standard specified does not have data in the standards directory'
        as_flux = data['FLUX'] * u.erg/ u.s / u.cm**2 / u.angstrom
        as_flux = as_flux.to(u.photon/u.s/u.angstrom/u.cm**2,
                             equivalencies=u.spectral_density(data['WAVELENGTH'] * u.AA))
        # Now calculate the approximate telescope area in cm**2
        telescope_area = np.pi * (8.1 * 100. / 2.)**2
        # Compute bandwidth per spectral pixel, assuming a fixed resolving power per
        # pixel of 200,000 and a total instrument throughput of 10%
        # FIXME: This is approximately 
        # applicable to GHOST, but could be more generally applicable and include
        # instrumental wavelength-dependent throughput.
        flux = as_flux * telescope_area * (data['WAVELENGTH'] / 2e5) * 0.1
        spectrum = np.array([data['WAVELENGTH'] / 1e4, np.maximum(flux,0)])
        return spectrum

        
    def simulate_observation(self, duration=0.0, output_prefix='test_',
                             spectrum_in=None, xshift=0.0, yshift=0.0,
                             rv_thar=0.0, flux=1, use_thar=True, res='high',
                             add_sky=True, radvel=0.0, thar_flatlamp=False,
                             flatlamp=False, obstype=None, objname=None,
                             seeing=0.8, sv_duration=10.0, sv_flux_profile=None,
                             data_label=1, binmode=(1, 1), utstart=None):
        """
        Simulate an observation with the whole instrument.
        This includes slit-viewing exposures, a blue and a red science exposure.

        Parameters
        ----------
        duration: float (optional)
            Duration of the exposure in seconds

        output_prefix: string (optional)
            Prefix for the output filename.

        spectrum_in: array of shape (2,n) where n is the number of
            spectral samples input spectrum.
            spectrum[0] is the wavelength array, spectrum[1] is the flux array,
            in units of recorded photons/spectral pixel/s, taking into account
            the full instrument at order centers (i.e. excluding the blaze
            function)

        xshift: float (optional)
            x-direction (along-slit) shift.

        yshift: float (optional)
            y-direction (spectral direction) shift.

        radvel: float (optional)
            Radial velocity in m/s for the target star with respect to
            the observer.

        rv_thar: float (optional)
            Radial velocity in m/s applied to the Thorium/Argon source.  It is
            unlikely that this is useful (use yshift instead for common shifts
            in the dispersion direction).

        flux: float (optional) Flux multiplier for the reference spectrum to
            give photons/pix/s.

        use_thar: bool (optional)
            Is the Thorium/Argon lamp in use?

        res: string (optional)
            Can be 'high' or 'std' for the resolution mode.

        add_sky: bool (optional)
            Is the sky background added to the frame?

        thar_flatlamp: bool (optional)
            Is the ThAr fiber illuminated with the flat lamp?

        flatlamp: bool (optional)
            Is the input spectrum evenly illuminating the slit, like a
            flat lamp would?

        obstype: string
            The value for the OBSTYPE FITS keyword

        obsname: string
            The name of the target to be created. Will populate the OBJECT
            keyword.

        seeing: float (optional)
            Determines the flux in each fiber (if flatlamp=False)

        sv_duration: the slit viewing camera exposure duration (in seconds)

        sv_flux_profile: array of shape (2, n) where n is the number of samples
            in the profile.  sv_flux_profile[0] is the sv exposure index array,
            and sv_flux_profile[1] gives the flux of that exposure, relative to
            the maximum possible (i.e. 0=no flux, 1=full flux).

        data_label: int (optional)
            Indicates which frame of a sequence this is; written to the DATALAB
            keyword
        """

        if utstart is None:
            utstart = datetime.datetime.utcnow()

        if (sv_duration > 0) and (duration > 0):
            self.slitv.set_exposure(
                sv_duration, duration, sv_flux_profile, utstart)
        else:
            self.slitv.set_exposure(0.0, 0, None, utstart)

        common_params = dict(
            duration=duration, output_prefix=output_prefix, gain=self.gain,
            spectrum_in=spectrum_in, xshift=xshift, yshift=yshift,
            radvel=radvel, rv_thar=rv_thar, flux=flux, rnoise=self.rnoise,
            bias_level=self.bias_level, overscan=self.overscan, add_sky=add_sky,
            namps=self.namps, use_thar=use_thar, res=res, utstart=utstart,
            thar_flatlamp=thar_flatlamp, seeing=seeing, data_label=data_label,
            flatlamp=flatlamp, obstype=obstype, objname=objname, binning=binmode)

        hdul_b = self.blue.simulate_frame(  # pylint: disable=star-args
            additive_noise=self.additive_noise['blue'],
            scaling=self.scaling['blue'], **common_params)
        hdul_r = self.red.simulate_frame(  # pylint: disable=star-args
            additive_noise=self.additive_noise['red'],
            scaling=self.scaling['red'], **common_params)
        hdul_s = self.slitv.save(output_prefix, obstype, res,
            data_label=data_label)  # noqa

        if self.split:
            return

        nextns_s = 0 if hdul_s is None else len(hdul_s)-1
        nextns_r = 0 if hdul_r is None else len(hdul_r)-1
        nextns_b = 0 if hdul_b is None else len(hdul_b)-1

        hdr = hdul_r[0].header
        del hdr['DETTYPE']
        del hdr['DETSIZE']
        del hdr['CAMERA']
        # hdr['EXTEND'] = True  # this has no effect for some reason
        hdr['NEXTEND'] = (nextns_r+nextns_b+nextns_s, 'Number of extensions')
        hdr['NREDEXP'] = (int(nextns_r/4), 'Number of red exposures')
        hdr['NBLUEEXP'] = (int(nextns_b/4), 'Number of blue exposures')
        hdr['NSLITEXP'] = (nextns_s, 'Number of slit-viewing exposures')
        hdulist = pf.HDUList([pf.PrimaryHDU(header=hdr)])
        if hdul_s is not None:
            hdulist += hdul_s[1:]
        hdulist += hdul_r[1:]
        hdulist += hdul_b[1:]

        print('Writing observation MEF to ' + output_prefix + 'MEF.fits')
        hdulist[0].header['EXTEND'] = True  # this works
        hdulist.writeto(output_prefix + 'MEF.fits', overwrite=True)
