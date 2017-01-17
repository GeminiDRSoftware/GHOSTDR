
"""This is a simple simulation code for GHOST or Veloce,
with a class ARM that simulates
a single arm of the instrument. The key default parameters
are hardwired for each named
configuration in the __init__ function of ARM.

Note that in this simulation code, the 'x' and 'y' directions are the
along-slit and dispersion directions respectively...
(similar to physical axes) but by convention, images are returned/displayed
 with a vertical slit and a horizontal dispersion direction.

For a simple simulation, run:

import pymfe

blue = pymfe.ghost.Arm('blue')

blue.simulate_frame()

TODO:
1) Add spectrograph aberrations (just focus and coma)
2) Add pupil illumination plus aberrations.
"""
import os
import numpy as np
from astrodata_GHOST.polyfit.polyspect import Polyspect


class Arm(Polyspect):
    """A class for each arm of the spectrograph. The initialisation
    function takes a series of strings representing the configuration.
    It can be "red" or "blue" for the arm (first string),
    and "std" or "high" for the mode (second string).
    """

    def __init__(self, arm='blue', mode='std'):
        """Initialisation function that sets all the mode specific parameters
        related to each configuration of the spectrograph.
        """
        # A lot of these parameters are yet unused.
        # These are a legacy of the original simulator and are left here
        # because they may become useful in the future.
        self.spect = 'ghost'
        self.arm = arm
        self.distance = 1000 / 52.67  # Distance in microns
        self.theta = 65.0  # Blaze angle
        self.assym = 1.0 / 0.41  # Magnification
        self.gamma = 0.56  # Echelle gamma
        self.nwave = 1e2  # Wavelengths per order for interpolation.
        self.f_col = 1750.6  # Collimator focal length.
        self.lenslet_high_size = 118.0  # Lenslet flat-to-flat in microns
        self.lenslet_std_size = 197.0  # Lenslet flat-to-flat in microns
        # When simulating the slit image, use this many microns per pixel
        self.microns_pix = 2.0
        # Number of microns in the slit image plane per arcsec
        self.microns_arcsec = 400.0
        self.im_slit_sz = 2048  # Size of the image slit size in pixels.
        # True if the spectral dispersion dimention is over the x (column) axis
        self.transpose = True
        self.mode = mode
        # This is the location of the model parameter files.
        # At this time this is not a brilliant idea as it is dependent on the
        # date of the correct version of the model and is a relative location
        # with respect to whichever module is imported. This relies on the
        # data being in a correct relative directory to polyfit
        # (or astrodata_GHOST, depending on whether the user imports this
        # directly).
        # Ideally we will want some sort of function that decides which model
        # to use for each data set, and overwrite this variable.
        # For now, it defaults to the latest one set manually.
        self.model_location = os.path.abspath(os.path.join(
            os.path.dirname(__file__),
            '../ADCONFIG_GHOST/lookups/GHOST/models/' +
            self.arm + '/161120/' + self.mode))
        if arm == 'red':
            # Additional slit rotation across an order needed to match Zemax.
            self.extra_rot = 3.0
            self.szx = 6144
            self.szy = 6160
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.drot = -2.0  # Detector rotation
            self.d_x = 1000 / 565.  # VPH line spacing
            self.theta_i = 30.0  # Prism incidence angle
            self.alpha1 = 0.0  # First prism apex angle
            self.alpha2 = 0.0  # Second prism apex angle
            self.m_min = 34
            self.m_max = 67
            self.m_ref = 50  # Reference order
            # Now put in the default fiber profile parameters for each mode.
            # These are used by the convolution function on polyspect
            # These were determined based on visual correspondence with
            # simulated data and may need to be revised once we have real
            # data. The same applies to the blue arm parameters.
            if self.mode == 'std':
                self.fiber_separation = 4.15
                self.profile_sigma = 1.1
            elif self.mode == 'high':
                self.fiber_separation = 2.49
                self.profile_sigma = 0.7
        elif arm == 'blue':
            # Additional slit rotation accross an order needed to match Zemax.
            self.extra_rot = 2.0
            self.szx = 4096
            self.szy = 4112
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.d_x = 1000 / 1137.  # VPH line spacing
            self.theta_i = 30.0  # Prism incidence angle
            self.drot = -2.0  # Detector rotation.
            self.alpha1 = 0.0  # First prism apex angle
            self.alpha2 = 0.0  # Second prism apex angle
            self.m_min = 63
            self.m_max = 95
            self.m_ref = 80  # Reference order
            # Now put in the default fiber profile parameters for each mode.
            # These are used by the convolution function on polyspect
            if self.mode == 'std':
                self.fiber_separation = 3.97
                self.profile_sigma = 1.1
            elif self.mode == 'high':
                self.fiber_separation = 2.53
                self.profile_sigma = 0.7
        else:
            print("Unknown spectrograph arm!")
            raise UserWarning
        # Now we determine the number of fibers based on mode.
        if mode == 'high':
            self.lenslet_width = self.lenslet_high_size
            self.nlenslets = 28
            # Set default profiles - object, sky and reference
            fluxes = np.zeros((self.nlenslets, 3))
            fluxes[2:21, 0] = 0.37
            fluxes[8:15, 0] = 0.78
            fluxes[11, 0] = 1.0
            # NB if on the following line, fluxes[2:,1]=1.0 is set, sky will be
            # subtracted automatically.
            fluxes[2 + 19:, 1] = 1.0
            fluxes[0, 2] = 1.0
        elif mode == 'std':
            self.lenslet_width = self.lenslet_std_size
            self.nlenslets = 17
            # Set default profiles - object 1, sky and object 2
            fluxes = np.zeros((self.nlenslets, 3))
            fluxes[0:7, 0] = 1.0
            fluxes[7:10, 1] = 1.0
            fluxes[10:, 2] = 1.0
        else:
            print("Unknown mode!")
            raise UserWarning

        Polyspect.__init__(self, self.mode, self.m_ref, self.szx, self.szy,
                           self.m_min, self.m_max, self.transpose,
                           self.extra_rot,self.nlenslets, self.fiber_separation,
                           self.profile_sigma)

    def make_lenslets(self, fluxes=[], seeing=0.8, llet_offset=0):
        """Make an image of the lenslets with sub-pixel sampling.
           This is potentially to become part of the simulator.

        Parameters
        ----------
        fluxes: float array (optional)
            Flux in each lenslet

        mode: string (optional)
            'high' or 'std', i.e. the resolving power mode of the spectrograph.
            Either mode or fluxes must be set.

        seeing: float (optional)
            If fluxes is not given, then the flux in each lenslet is defined
            by the seeing.

        llet_offset: int
            Offset in lenslets to apply to the input spectrum"""
        print("Computing a simulated slit image...")
        szx = self.im_slit_sz
        szy = 256
        fillfact = 0.98
        s32 = np.sqrt(3) / 2
        hex_scale = 1.15
        # equivalent to a 1 degree FWHM for an f/3 input ???
        conv_fwhm = 30.0

        if self.mode == 'std':
            yoffset = (self.lenslet_width / self.microns_pix / hex_scale *
                       np.array([0, -s32, s32, 0, -s32, s32, 0])).astype(int)
            xoffset = (self.lenslet_width / self.microns_pix /
                       hex_scale * np.array([-1, -0.5, -0.5, 0,
                                             0.5, 0.5, 1.0])).astype(int)
        elif self.mode == 'high':
            yoffset = (self.lenslet_width / self.microns_pix /
                       hex_scale * s32 *
                       np.array([-2, 2, -2, -1, -1, 0, -1, -1, 0, 0, 0,
                                 1, 1, 0, 1, 1, 2, -2, 2])).astype(int)
            xoffset = (self.lenslet_width / self.microns_pix / hex_scale *
                       0.5 * np.array([-2, 0, 2, -3, 3, -4, -1, 1, -2, 0, 2,
                                       -1, 1, 4, -3, 3, -2, 0, 2])).astype(int)
        else:
            print("Error: mode must be standard or high")

        # Some preliminaries...
        cutout_hw = int(self.lenslet_width / self.microns_pix * 1.5)
        im_slit = np.zeros((szy, szx))
        x_values = np.arange(szx) - szx / 2.0
        y_values = np.arange(szy) - szy / 2.0
        xygrid = np.meshgrid(x_values, y_values)
        # r and wr enable the radius from the lenslet center to be indexed
        radius = np.sqrt(xygrid[0]**2 + xygrid[1]**2)
        wradius = np.where(radius < 2 * self.lenslet_width / self.microns_pix)
        # g is a Gaussian used for FRD
        gauss = np.exp(-radius**2 / 2.0 / (conv_fwhm / self.microns_pix / 2.35)**2)
        gauss = np.fft.fftshift(gauss)
        gauss /= np.sum(gauss)
        gft = np.conj(np.fft.rfft2(gauss))
        pix_size_slit = self.px_sz * \
            (self.f_col / self.assym) / self.f_cam * 1000.0 / self.microns_pix
        pix = np.zeros((szy, szx))
        pix[np.where((np.abs(xygrid[0]) < pix_size_slit / 2) *
                     (np.abs(xygrid[1]) < pix_size_slit / 2))] = 1
        pix = np.fft.fftshift(pix)
        pix /= np.sum(pix)
        pix_ft = np.conj(np.fft.rfft2(pix))
        # Create some hexagons. We go via a "cutout" for efficiency.
        h_cutout = self.hexagon(
            szy, self.lenslet_width / self.microns_pix * fillfact / hex_scale)
        hbig_cutout = self.hexagon(
            szy, self.lenslet_width / self.microns_pix * fillfact)
        height = np.zeros((szy, szx))
        hbig = np.zeros((szy, szx))
        height[:, szx / 2 - szy / 2:szx / 2 + szy / 2] = h_cutout
        hbig[:, szx / 2 - szy / 2:szx / 2 + szy / 2] = hbig_cutout
        if len(fluxes) != 0:
            # If we're not simulating seeing, the image-plane is uniform, and
            # we only use the values of "fluxes" to scale the lenslet fluxes.
            image = np.ones((szy, szx))
            # Set the offsets to zero because we may be simulating a single
            # Th/Ar lenslet and not starlight (from the default xoffset etc)
            xoffset = np.zeros(len(fluxes), dtype=int)
            yoffset = np.zeros(len(fluxes), dtype=int)
        else:
            # If we're simulating seeing, create a Moffat function as our input
            # profile,but just make the lenslet fluxes uniform.
            image = np.zeros((szy, szx))
            im_cutout = self.moffat2d(
                szy, seeing * self.microns_arcsec / self.microns_pix / 2)
            image[:, szx / 2 - szy / 2:szx / 2 + szy / 2] = im_cutout
            fluxes = np.ones(len(xoffset))

        # Go through the flux vector and fill in each lenslet.
        for i in range(len(fluxes)):
            im_one = np.zeros((szy, szx))
            im_cutout = np.roll(
                np.roll(image, yoffset[i], axis=0), xoffset[i], axis=1) * height
            im_cutout = im_cutout[szy / 2 - cutout_hw:szy / 2 +
                                  cutout_hw, szx / 2 - cutout_hw:szx /
                                  2 + cutout_hw]
            prof = self.azimuthal_average(
                im_cutout, returnradii=True, binsize=1)
            prof = (prof[0], prof[1] * fluxes[i])
            xprof = np.append(np.append(0, prof[0]), np.max(prof[0]) * 2)
            yprof = np.append(np.append(prof[1][0], prof[1]), 0)
            im_one[wradius] = np.interp(radius[wradius], xprof, yprof)
            im_one = np.fft.irfft2(np.fft.rfft2(im_one) * gft) * hbig
            im_one = np.fft.irfft2(np.fft.rfft2(im_one) * pix_ft)
            # !!! The line below could add tilt offsets...
            # important for PRV simulation !!!
            # im_one = np.roll(np.roll(im_one, tilt_offsets[0,i], axis=1),
            # tilt_offsets[1,i], axis=0)*hbig
            the_shift = int((llet_offset + i - self.nlenslets / 2.0) *
                            self.lenslet_width / self.microns_pix)
            im_slit += np.roll(im_one, the_shift, axis=1)
        return im_slit

    #-- This is a collection of functions that were part of optics.py --
    @staticmethod
    def hexagon(dim, width):
        """This function creates a hexagon. Originally from opticstools.

        Parameters
        ----------
        dim: int
            Size of the 2D array
        width: int
            flat-to-flat width of the hexagon

        Returns
        -------
        pupil: float array (sz,sz)
            2D array hexagonal pupil mask
        """
        x_values = np.arange(dim) - dim / 2.0
        xygrid = np.meshgrid(x_values, x_values)
        xxg = xygrid[1]
        yyg = xygrid[0]
        w_condition = np.where((yyg < width / 2) * (yyg > (-width / 2)) *
                               (yyg < (width - np.sqrt(3) * xxg)) *
                               (yyg > (-width + np.sqrt(3) * xxg)) *
                               (yyg < (width + np.sqrt(3) * xxg)) *
                               (yyg > (-width - np.sqrt(3) * xxg)))
        hexgon = np.zeros((dim, dim))
        hexgon[w_condition] = 1.0
        return hexgon

    def moffat(theta, halfwidth, beta=4.0):
        """This creates a moffatt function for simulating seeing.
        The output is an array with the same dimensions as theta.
        Total Flux" is set to 1 - this only applies if sampling
        of thetat is 1 per unit area (e.g. arange(100)).

        From Racine (1996), beta=4 is a good approximation for seeing

        Parameters
        ----------
        theta: float or float array
            Angle at which to calculate the moffat profile (same units as hw)
        hw: float
            Half-width of the profile
        beta: float
            beta parameters

        """
        denom = (1 + (2**(1.0 / beta) - 1) * (theta / halfwidth)**2)**beta
        return (2.0**(1.0 / beta) - 1) * (beta - 1) / np.pi / halfwidth**2 / denom

    @staticmethod
    def moffat2d(szbase, halfwidth, beta=4.0):
        """A 2D version of a moffat function
        """
        x_values = np.arange(szbase) - szbase / 2.0
        xygrid = np.meshgrid(x_values, x_values)
        radius = np.sqrt(xygrid[0]**2 + xygrid[1]**2)
        mof = moffat(radius, halfwidth, beta=beta)
        return mof

    @staticmethod
    def azimuthal_average(image, center=None, stddev=False, returnradii=False,
                          return_nr=False, binsize=0.5, weights=None,
                          steps=False, interpnan=False, left=None,
                          right=None, return_max=False):
        """
        Calculate the azimuthally averaged radial profile.
        NB: This was found online and should be properly credited!
        Modified by MJI

        image - The 2D image
        center - The [x,y] pixel coordinates used as the center.
                 The default is None, which then uses the center of the image
                 (including fractional pixels).
        stddev - if specified, return the azimuthal standard deviation instead
                 of the average
        returnradii - if specified, return (radii_array,radial_profile)
        return_nr   - if specified, return number of pixels per radius *and*
                      radius
        binsize - size of the averaging bin.  Can lead to strange results if
            non-binsize factors are used to specify the center and the binsize
            is too large
        weights - can do a weighted average instead of a simple average if
                  this keyword parameter is set.
                  weights.shape must = image.shape.
                  weighted stddev is undefined, so don't set weights and stddev.
        steps - if specified, will return a double-length bin array and radial
            profile so you can plot a step-form radial profile (which more
            accurately represents what's going on)
        interpnan - Interpolate over NAN values, i.e. bins where there is no
                    data?
            left,right - passed to interpnan; they set the extrapolated values
        return_max - (MJI) Return the maximum index.

        If a bin contains NO DATA, it will have a NAN value because of the
        divide-by-sum-of-weights component.  I think this is a useful way to
        denote lack of data, but users let me know if an alternative is
        prefered...

        """
        # Calculate the indices from the image
        y, x = np.indices(image.shape)

        if center is None:
            center = np.array([(x.max() - x.min()) / 2.0,
                               (y.max() - y.min()) / 2.0])

        r = np.hypot(x - center[0], y - center[1])

        if weights is None:
            weights = np.ones(image.shape)
        elif stddev:
            raise ValueError("Weighted standard deviation is not defined.")

        # the 'bins' as initially defined are lower/upper bounds for each bin
        # so that values will be in [lower,upper)
        nbins = int(np.round(r.max() / binsize) + 1)
        maxbin = nbins * binsize
        bins = np.linspace(0, maxbin, nbins + 1)
        # but we're probably more interested in the bin centers than their left
        # or right sides...
        bin_centers = (bins[1:] + bins[:-1]) / 2.0

        # Find out which radial bin each point in the map belongs to
        whichbin = np.digitize(r.flat, bins)

        # how many per bin (i.e., histogram)?
        # there are never any in bin 0, because the lowest index returned by
        # digitize is 1
        nr = np.bincount(whichbin)[1:]

        # recall that bins are from 1 to nbins (which is expressed in array
        # terms by arange(nbins)+1 or xrange(1,nbins+1) )
        # radial_prof.shape = bin_centers.shape

        if stddev:
            radial_prof = np.array([image.flat[whichbin == b].std()
                                    for b in xrange(1, nbins + 1)])
        elif return_max:
            radial_prof = np.array([np.append(
                (image * weights).flat[whichbin == b], -np.inf).max()
                                    for b in xrange(1, nbins + 1)])
        else:
            radial_prof = np.array([(image * weights).flat[whichbin == b].sum()
                                    / weights.flat[whichbin == b].sum()
                                    for b in xrange(1, nbins + 1)])

        #import pdb; pdb.set_trace()

        if interpnan:
            radial_prof = np.interp(bin_centers,
                                    bin_centers[radial_prof == radial_prof],
                                    radial_prof[radial_prof == radial_prof],
                                    left=left, right=right)

        if steps:
            xarr = np.array(zip(bins[:-1], bins[1:])).ravel()
            yarr = np.array(zip(radial_prof, radial_prof)).ravel()
            return xarr, yarr
        elif returnradii:
            return bin_centers, radial_prof
        elif return_nr:
            return nr, bin_centers, radial_prof
        else:
            return radial_prof
