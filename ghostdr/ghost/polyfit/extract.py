"""
Given (x,wave,matrices, slit_profile), extract the flux from each order.
"""

from __future__ import division, print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
import matplotlib.cm as cm
import warnings
import scipy.ndimage as ndimage
from scipy import optimize
from scipy.interpolate import CubicSpline
from gempy.library import astrotools, matching
from geminidr.gemini.lookups import DQ_definitions as DQ


def find_additional_crs(phi, col_data, col_inv_var, snoise=0.1, sigma=6,
                        noise_model=None, reject=True, debug=False, pixel=None):
    """
    Utility function to search for additional cosmic rays.

    Pixels are considered to be cosmic rays if their pixel value exceeds

    .. math::
        \\hat{y} + n_\\sigma \\times \\sqrt{var}\\textrm{,}

    where:

    - :math:`n_\\sigma` is the number of standard deviations (``nsigma``);
    - :math:`\\hat{y}` is the maximum value of a linear regression model
      fitted to the column data;
    - :math:`var` is an effective variance for the column data, computed as:

      .. math::
          var = \\frac{1}{\\verb+col_inv_var+ + ( \\verb+snoise+ \\times \\verb+coldata+ )^2}
    
    Parameters
    ----------
    phi: :obj:`numpy.ndarray`
        A (nobj :math:`\\times` npix) model PSF array.
    col_data: :obj:`numpy.ndarray`
        A (npix) data array
    col_inv_var: :obj:`numpy.ndarray`
        A (npix) data variance array
    s_noise : float
        A noise factor to be used in the cosmic ray detection calculation.
    sigma : int
        Number of standard deviations permitted before a pixel is flagged as
        bad.
    noise_model: astropy.modeling.Model object
        a model to predict the variance of a pixel with a given signal
    reject: bool
        actually reject pixels? (rather than just fit spatial profiles)
    debug : bool, default False
        If True, show plots
    pixel : tuple
        order and location of this pixel (for additional debugging)
        
    Returns
    -------
    new_bad: :obj:`numpy.ndarray`
        Array of indices of bad pixels in ``col_data``.
    """
    spat_conv_weights = np.array([0.25, .5, .25])
    # Don't try to fit points which aren't expected to have any signal
    orig_good = np.logical_and(col_inv_var > 0,
                               phi.sum(axis=0) > 0)

    # If we don't have enough pixels to fit we can't reject
    if orig_good.sum() < phi.shape[0]:
        result = optimize.lsq_linear(phi.T, col_data, bounds=(0, np.inf))
        return [], result.x

    good = orig_good.copy()
    if debug:
        print("\n")
        print("OBJECT PROFILES (phi)")
        print(phi)
    # Don't use the weights initially because these are smoothed
    # and small CRs get missed
    try:
        result = optimize.lsq_linear(phi.T[good], col_data[good], bounds=(0, np.inf))
    except:
        print(f"INITIAL FITTING FAILURE AT PIXEL {pixel}")
        print(col_data[good])
        for p in phi:
            print([p[good]])
        raise
    if not reject:  # we have a result, initial VAR array OK for weights
        return [], result.x

    model = np.dot(result.x, phi)
    inv_var_use = 1. / np.maximum(noise_model(model) + (snoise * model) ** 2, 0)
    nit = 0
    while good.sum() >= phi.shape[0]:
        # CJS problem with using weights and rejecting low pixels
        # Further investigation needed... can we be sure there are
        # no low pixels?
        A = phi * np.sqrt(inv_var_use)
        b = col_data * np.sqrt(inv_var_use)
        try:
            result = optimize.lsq_linear(A.T[good], b[good], bounds=(0, np.inf))
        except:
            print(f"FITTING FAILURE AT PIXEL {pixel} ITERATION {nit}")
            print(col_data[good])
            for p in phi:
                print(p[good])
            print("-"*60)
            print(inv_var_use)
            print(model)
            print("-"*60)
            print(b[good])
            for p in A:
                print(p[good])
            raise

        model = np.dot(result.x, phi)
        if debug:
            print("RAW COLUMN DATA")
            print(col_data)
            print("RESULT", result.x)
            print(model)
            print(good)
            print("-" * 60)
        var_use = ndimage.convolve1d(np.maximum(noise_model(model) + (snoise * model) ** 2, 0), spat_conv_weights)
        inv_var_use = np.where(var_use > 0, 1. / var_use, 0)
        deviation = (col_data - model) * np.sqrt(inv_var_use)
        if debug:
            print(np.sqrt(var_use))
            print(deviation)
            print("-" * 60)
        worst_offender = np.argmax(abs(deviation)[good])
        if abs(deviation[good][worst_offender]) > sigma:
            good[np.where(good)[0][worst_offender]] = False
        #new_bad = (abs(col_data - model) > max_deviation) & good
        #if new_bad.any():
        #    good &= ~new_bad
        else:
            break
        nit += 1

    new_bad = np.where(np.logical_and.reduce(
        [col_inv_var > 0, abs(deviation) > sigma, orig_good]))[0]
    #print("RESULT", *pixel, result.x)

    # CJS 20230120: this allows a little extra leeway for vertical CCD bleed
    #limit = ndimage.maximum_filter(y_hat + nsigma * np.sqrt(var_use),
    #                               size=3, mode='constant')
    #limit = y_hat + nsigma * np.sqrt(var_use)
    #new_bad = np.where(col_data > limit)[0]
    #if debug and len(new_bad) > 0:
    if debug:
        print("BAD", new_bad, phi.shape)
        # Sky is always last profile in phi?
        summed_fluxes = [np.sum(col_data[p>0]) for p in phi[:-1]]
        print("SUMMED FLUXES IN OBJECT FIBRES (not sky subtracted)",
              summed_fluxes)
        plt.ioff()
        plt.plot(col_data, 'k-', label='data')
        # max allowed value
        plt.plot(model + sigma * np.sqrt(var_use), 'r-', label='limit')
        plt.plot(model, 'b-', label='model')
        for ppp, amp in zip(phi, result.x):
            plt.plot(ppp * amp, ':')
        plt.show()
        plt.ion()

    return new_bad, result.x


def subtract_scattered_light(data, mask):
    """
    Remove (approximately) scattered light on the detector.

    .. note::
        This implementation is a no-op function, which will simply return
        ``data``.

    Subtract a simple linear approximation to the scattered light
    on the detector. This is a place-holder function, which can be upgraded to 
    higher order scattered light subtraction, depending on the performance of 
    the spectrograph.
    
    Parameters
    ----------
    data: :obj:`numpy.ndarray`
        The data to be corrected for scattered light.
    
    mask: :obj:`numpy.ndarray`
        A mask which is True wherever there is only scattered light and no data.
    
    Returns
    -------
    data: :obj:`numpy.ndarray`
        A scattered-light corrected image.
    
    """
    return data


class Extractor(object):
    """
    A class to extract data for each arm of the spectrograph.

    The extraction is defined by 3 key parameters: an ``x_map``, which is
    equivalent to 2dFDR's tramlines and contains a physical x-coordinate for
    every y (dispersion direction) coordinate and order, and a ``w_map``,
    which is the wavelength corresponding to every y (dispersion direction)
    coordinate and order.

    The slit parameters are defined by the slitview_instance. If a different
    slit profile is to be assumed, then another (e.g. simplified or fake)
    slitview instance should be made, rather than modifying this class.

    Attributes
    ----------
    arm: :any:`polyspect.Polyspect`
        This defines, e.g., whether the camera is ``"red"`` or ``"blue"``

    slitview: :any:`slitview.SlitView`
        This defines, e.g., whether the mode is ``"std"`` or ``"high"``

    gain: float, optional
        gain in electrons per ADU. From fits header.

    rnoise: float, optional
        Expected readout noise in electrons. From fits header.

    badpixmask: :obj:`numpy.ndarray`, optional
        A data quality plane, which evaluates to False (i.e. 0) for good
        pixels.

    vararray: :obj:`numpy.ndarray`, optional
        A variance array.

    transpose: bool , optional
        Do we transpose the data before extraction? Default is ``False``.

    cr_flag: integer, optional
        When we flag additional cosmic rays in the badpixmask, what value
        should we use? Default is ``8``.
    """
    def __init__(self, polyspect_instance, slitview_instance,
                 gain=1.0, rnoise=3.0, cr_flag=DQ.cosmic_ray,
                 badpixmask=None, transpose=False,
                 vararray=None):
        self.arm = polyspect_instance
        self.slitview = slitview_instance
        self.transpose = transpose
        self.gain = gain
        self.rnoise = rnoise
        self.vararray = vararray
        self.badpixmask = badpixmask
        self.cr_flag = cr_flag

        # FIXME: This warning could probably be neater.
        if not isinstance(self.arm.x_map, np.ndarray):
            raise UserWarning('Input polyspect_instance requires'
                              'spectral_format_with matrix to be run.')

        # To aid in 2D extraction, let's explicitly compute the y offsets
        # corresponding to these x offsets...
        # The "matrices" map pixels back to slit co-ordinates.
        ny = self.arm.x_map.shape[1]
        nm = self.arm.x_map.shape[0]
        self.slit_tilt = np.zeros((nm, ny))
        for i in range(nm):
            for j in range(ny):
                invmat = np.linalg.inv(self.arm.matrices[i, j])
                # What happens to the +x direction?
                x_dir_map = np.dot(invmat, [1, 0])
                self.slit_tilt[i, j] = x_dir_map[1] / x_dir_map[0]

    def bin_models(self):
        """
        Bin the models to match data binning.

        Function used to artificially bin the models so that they apply to
        whatever binning mode the data are. This requires knowledge of the 
        x and y binning from the arm class, which is assumed this class
        inherits. 
        
        The binning is done as a running average, in which the
        values for each binned pixel are assumed to be equivalent to the average
        value of all physical pixels that are part of the binned pixel.

        Returns
        -------
        x_map: :obj:`numpy.ndarray`
            Binned version of the x map model
        wmap: :obj:`numpy.ndarray`
            Binned version of the wavelength map model
        blaze: :obj:`numpy.ndarray`
            Binned version of the blaze map model
        matrices: :obj:`numpy.ndarray`
            Binned version of the matrices array
        """
        if self.arm.xbin == 1 and self.arm.ybin == 1:
            return self.arm.x_map, self.arm.w_map, self.arm.blaze, \
                   self.arm.matrices
        # Start by getting the order number. This should never change.
        n_orders = self.arm.x_map.shape[0]

        x_map = self.arm.x_map.copy()
        w_map = self.arm.w_map.copy()
        blaze = self.arm.blaze.copy()
        matrices = self.arm.matrices.copy()
        # The best way to do this is to firstly do all the ybinning, and then do
        # the x binning
        if self.arm.ybin > 1:
            # Now bin the x_map, firstly in the spectral direction
            # We do this by reshaping the array by adding another dimension of
            # length ybin and then averaging over this axis
            x_map = np.mean(x_map.reshape(n_orders,
                                          int(self.arm.szy / self.arm.ybin),
                                          self.arm.ybin), axis=2)

            # Now do the same for the wavelength scale and blaze where necessary
            w_map = np.mean(w_map.reshape(n_orders,
                                          int(self.arm.szy / self.arm.ybin),
                                          self.arm.ybin), axis=2)

            blaze = np.mean(blaze.reshape(n_orders,
                                          int(self.arm.szy / self.arm.ybin),
                                          self.arm.ybin), axis=2)
            # The matrices are a bit harder to work with, but still the same
            # principle applies.
            matrices = np.mean(matrices.reshape(n_orders,
                                                int(
                                                    self.arm.szy /
                                                    self.arm.ybin),
                                                self.arm.ybin, 2, 2), axis=2)

        if self.arm.xbin > 1:
            # Now, naturally, the actualy x values must change according to the
            # xbin
            #x_map /= self.arm.xbin
            x_map = (x_map + 0.5) / self.arm.xbin - 0.5
            # The w_map and blaze remain unchanged by this. 

        # Now we must modify the values of the [0,0] and [1,1] elements of
        # each matrix according to the binning to reflect the now size of
        # binned pixels.
        rescale_mat = np.array([[self.arm.xbin, 0], [0, self.arm.ybin]])
        matrices = np.dot(matrices, rescale_mat)

        # Set up convenience local variables
        ny = x_map.shape[1]
        nm = x_map.shape[0]

        self.slit_tilt = np.zeros((nm, ny))
        for i in range(nm):
            for j in range(ny):
                invmat = np.linalg.inv(matrices[i, j])
                # What happens to the +x direction?
                x_dir_map = np.dot(invmat, [1, 0])
                self.slit_tilt[i, j] = x_dir_map[1] / x_dir_map[0]

        # Hack to see if this is the only problem
        #matrices *= 0.99

        return x_map, w_map, blaze, matrices

    def make_pixel_model(self, input_image=None):
        """
        Based on the xmod and the slit viewer image, create a complete model image, 
        where flux versus wavelength pixel is constant. As this is designed for 
        comparing to flats, normalisation is to the median of the non-zero pixels in the
        profile.
        
        Parameters
        ----------
        input_iage: :obj:`numpy.ndarray`, optional
            Image data, transposed so that dispersion is in the "y" direction.
            If this is given, then the pixel model is scaled according to the input flux
            for every order and wavelength. Note that this isn't designed to reproduce
            dual-object or object+sky data.
        
        Returns
        -------
        model: :obj:`numpy.ndarray`
            An image the same size as the detector.
        """
        # Adjust in case of y binning (never for flats, which is what this function
        # is primarily designed for.
        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            raise RuntimeError('Extraction failed, unable to bin models.')
        #print("XMAP")
        #for i, x in enumerate(x_map[12]):
        #    print(f"{i:4d} {x + self.arm.szx // 2}")

        # Key array constants
        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx / self.arm.xbin)
                             
        profiles = [self.slitview.slit_profile(arm=self.arm.arm)]
        
        n_slitpix = profiles[0].shape[0]
        profile_y_microns = (np.arange(n_slitpix) -
                             n_slitpix / 2 + 0.5) * self.slitview.microns_pix
        
        if self.transpose:
            pixel_model = np.zeros( (ny, nx) ) 
        else:
            pixel_model = np.zeros( (nx, ny) )
        
        # Loop through all orders then through all y pixels.        
        print("    Creating order ", end="")
        for i in range(nm):
            print(f"{self.arm.m_min+i}...", end="")
            sys.stdout.flush()

            # Create an empty model array. Base the size on the largest
            # slit magnification for this order.
            nx_cutout = int(np.ceil(self.slitview.slit_length / np.min(
                matrices[i, :, 0, 0])))

            for j in range(ny):
                # Create an array of y pixels based on the scale for this
                # order and pixel.

                # Create our column cutout for the data and the PSF, then put
                # it back into our 2D array. 
                
                # FIXME: Is "round" most correct on the next line???
                # FIXME: Profiles should be convolved by a detector pixel,
                # but binning has to be taken into account properly!
                #x_ix = int(np.round(
                #    x_map[i, j])) - nx_cutout // 2 + \
                #       np.arange(nx_cutout, dtype=int) + nx // 2

                # CRH 20220825:  Convolve the slit to detector pixels 
                # by interpolating the slit profile and integrating it 
                # across the detector pixels

                #n_slit_sample = int(np.round(1.0/(profile_y_pix[1] - \
                #    profile_y_pix[0])))*2
                #x_ix_sample = int(np.round(
                #    x_map[i, j])) - nx_cutout // 2 + \
                #       np.arange(nx_cutout*n_slit_sample,
                #                 dtype=int)/n_slit_sample + nx // 2

                # Remove? These were old interpolations of the slit to the
                # detector

                # Depending on the slit orientation, we have to make sure 
                # that np.interp is fed a monotonically increasing x vector.
                #if matrices[i, j, 0, 0] < 0:
                #    phi = np.interp(
                #        x_ix - x_map[i, j] - nx // 2, profile_y_pix[::-1],
                #        profile[::-1])
                #else:
                #    phi = np.interp(
                #        x_ix - x_map[i, j] - nx // 2, profile_y_pix,
                #        profile)

                #interp_prof = np.interp(
                #    x_ix_sample - x_map[i, j] - nx // 2, profile_y_pix,
                #    profile)

                #phi = interp_prof.reshape(x_ix.shape[0],n_slit_sample).sum(axis=1)

                # Normalise to the median of the non-zero pixels. This is neater
                # if it normalises to the median, but normalisation can occur 
                # later, and shouldn't occur wavelength by wavelength.
                #phi /= np.sum(phi[phi != 0])

                x_ix, phi, profiles = resample_slit_profiles_to_detector(
                    profiles, profile_y_microns, x_map[i, j] + nx // 2,
                    detpix_microns=matrices[i, j, 0, 0])

                # This isn't perfect as it gets the top pixel wrong if the
                # trace goes over the top but it's OK for this purpose
                if self.transpose:
                    pixel_model[j, np.minimum(x_ix, nx-1)] = phi[0]
                else:
                    pixel_model[np.minimum(x_ix, nx-1), j] = phi[0]
                    
        return pixel_model

    def new_extract(self, data=None, correct_for_sky=True, use_sky=True,
                    optimal=False, find_crs=True, snoise=0.1, sigma=6,
                    used_objects=[0, 1], debug_cr_pixel=None,
                    correction=None):
        from .extractum import Extractum

        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            raise RuntimeError('Extraction failed, unable to bin models.')

        # Set up convenience local variables
        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx / self.arm.xbin)

        # Our profiles... we re-extract these in order to include the centroids
        # Centroids is the offset in pixels along the short axis of the pseudoslit
        profile, centroids = self.slitview.slit_profile(arm=self.arm.arm,
                                                        return_centroid=True)

        # Allow us to compute flux-weighted transverse position
        slitview_profiles = [profile, profile*centroids]

        def transverse_positions(slitview_profiles, profile_center=None,
                                 detpix_microns=None):
            """Returns offset along short axis of pseudoslit in microns"""
            x_ix, phi, _ = resample_slit_profiles_to_detector(
                slitview_profiles, profile_y_microns=profile_y_microns,
                profile_center=profile_center, detpix_microns=detpix_microns)
            return x_ix, astrotools.divide0(phi[1] / phi[0]) * self.slitview.microns_pix

        profiles = self.slitview.object_slit_profiles(
            arm=self.arm.arm, correct_for_sky=correct_for_sky,
            used_objects=used_objects, append_sky=use_sky or find_crs
        )

        # Number of "objects" and "slit pixels"
        no = profiles.shape[0]
        if find_crs and not use_sky:
            no -= 1
        n_slitpix = profiles.shape[1]
        profile_y_microns = (np.arange(n_slitpix) -
                             n_slitpix / 2 + 0.5) * self.slitview.microns_pix

        m_init = models.Polynomial1D(degree=1)
        fit_it = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip)
        good = ~np.isinf(self.vararray)
        if self.badpixmask is not None:
            good &= (self.badpixmask & DQ.not_signal) == 0
        noise_model, _ = fit_it(m_init, data[good].ravel(), self.vararray[good].ravel())
        if noise_model.c0 <= 0:
            print("Problem with read noise estimate!")
            print(noise_model)
            noise_model.c0 = 4  # just something for now

        pixel_inv_var = 1. / self.vararray

        #### STUFF ABOUT CONVOLUTION HERE ####

        # Loop through all orders then through all y pixels.
        if find_crs:
            print("    Extracting order ", end="")
            for i in range(nm):
                print(f"{self.arm.m_min+i}...", end="")
                sys.stdout.flush()

                for j in range(ny):
                    debug_this_pixel = debug_cr_pixel in [(self.arm.m_min+i, j)]

                    x_ix, phi, profiles = resample_slit_profiles_to_detector(
                        profiles, profile_y_microns, x_map[i, j] + nx // 2,
                        detpix_microns=matrices[i, j, 0, 0], debug=debug_this_pixel)

                    # Deal with edge effects...
                    ww = np.logical_or(x_ix >= nx, x_ix < 0)
                    x_ix[ww] = 0
                    phi /= phi.sum(axis=1)[:, np.newaxis]

                    _slice = (j, x_ix) if self.transpose else (x_ix, j)
                    col_data = data[_slice]
                    col_inv_var = pixel_inv_var[_slice]
                    badpix = (np.zeros_like(col_data, dtype=bool)
                              if self.badpixmask is None else self.badpixmask[_slice].astype(bool))
                    badpix[ww] = True
                    xtr = Extractum(phi, col_data, col_inv_var, badpix,
                                    noise_model=noise_model, pixel=(self.arm.m_min+i, j))


    def one_d_extract(self, data=None, fl=None, correct_for_sky=True,
                      use_sky=True, optimal=False, find_crs=True,
                      snoise=0.1, sigma=6,
                      used_objects=[0,1], vararray=None,
                      debug_cr_pixel=None, correction=None):
        """
        Extract flux by integrating down columns.

        This function extracts flux information by integrating down columns
        (the ``y`` direction), using an optimal extraction method.

        .. note::

            Given that some of this code is in common with two_d_extract, the
            routines could easily be merged; however, that would make
            one_d_extract less readable.

        Parameters
        ----------
        data: :obj:`numpy.ndarray` , optional
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        fl: string , optional
            A fits file with conventional row/column directions containing the
            data to be extracted.

        correct_for_sky: bool, optional
            Do we correct the object slit profiles for sky? Should be yes for
            objects and no for flats/arcs.
            
        use_sky: book, optional
            In extraction, do we use sky (and therefore self-subtract)?

        optimal: bool
            perform optimal (rather than uniform) extraction?

        find_crs : bool
            whether or not to call :any:`find_additional_crs`

        debug_cr_pixel : 2-tuple
            (order, pixel) of extraction to plot for debugging purposes

        snoise : float
            linear fraction of signal to add to noise estimate for CR flagging

        sigma : float
            number of standard deviations for identifying discrepant pixels

        vararray : :obj:`numpy.ndarray` , optional
            If given, the instance's `vararray` attribute will be updated
            to hold this array.

        Raises
        ------
        ValueError
            If neither ``data`` nor ``fl`` are passed.

        Returns
        -------
        extracted_flux: :obj:`numpy.ndarray`
            Extracted fluxes as a function of pixel along the spectral direction
        extracted_var: :obj:`numpy.ndarray`
            Extracted variance as a function of pixel along the spectral
            direction
        extraction_weights: :obj:`numpy.ndarray`
            Extraction weights as a function of pixel along the spectral
            direction

        """
        #Lets keep track of the number of additional cosmic rays as an object property
        self.num_additional_crs = 0

        # Update the instance .vararray, if a new one is passed
        if vararray is not None:
            self.vararray = vararray
        
        if data is None:
            if fl is None:
                raise ValueError("Must input data or file")
            else:
                data = pyfits.getdata(fl)

        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            raise RuntimeError('Extraction failed, unable to bin models.')

        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx / self.arm.xbin)

        # Our profiles...
        # FIXME: Consider carefully whether there is a way to extract x-centroids
        # as well for PRV, as part of slitim_offsets below.
        profiles = self.slitview.object_slit_profiles(
            arm=self.arm.arm, correct_for_sky=correct_for_sky,
            used_objects=used_objects, append_sky=use_sky or find_crs
        )

        # Number of "objects" and "slit pixels"
        no = profiles.shape[0]
        if find_crs and not use_sky:
            no -= 1
        n_slitpix = profiles.shape[1]
        profile_y_microns = (np.arange(n_slitpix) -
                             n_slitpix / 2 + 0.5) * self.slitview.microns_pix

        # To consider slit tilt for a 1D extraction, we need to know the profile
        # centroids in the "y" direction. i.e. In principle, we could modify the
        # wavelength scale for each object based on this. If 2D extraction works
        # well, such an approach is not needed, but lets keep the idea of this
        # code here for now. (code deleted but comment left as reminder)

        # Our extracted arrays, and the weights array
        extracted_flux = np.zeros((nm, ny, no), dtype=np.float32)
        extracted_var = np.zeros((nm, ny, no), dtype=np.float32)
        extraction_weights = np.zeros((no, nx, ny), dtype=np.float32)

        m_init = models.Polynomial1D(degree=1)
        fit_it = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip)
        good = ~np.isinf(vararray)
        if self.badpixmask is not None:
            good &= (self.badpixmask & DQ.not_signal) == 0
        noise_model, _ = fit_it(m_init, data[good].ravel(), vararray[good].ravel())
        if noise_model.c0 <= 0:
            print("Problem with read noise estimate!")
            print(noise_model)
            noise_model.c0 = 4  # just something for now

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel inverse variance.
        # This really should come from an input "vararray" because of differing
        # gain and readout noise parameters for different amplifiers, and known
        # bad pixels.
        if self.vararray is None:
            pixel_inv_var = 1.0 / (np.maximum(data, 0) /
                                   self.gain + self.rnoise ** 2)
        else:
            pixel_inv_var = 1.0 / self.vararray

        # Now, smooth filter this variance plane for the purposes of not allowing tilted
        # slits or line profiles to affect the extraction. Simply convolve the inverse
        # variance by a profile 7 wide in the spectral direction and preserve bad pixels. 
        # This ensures that line profile as extracted will not be SNR-dependent
        # (optimal extraction is only independent of line profiles for Gaussian 2D 
        # PSFs)
        spec_conv_weights = np.hanning(15)[1:-1]
        spec_conv_weights /= np.sum(spec_conv_weights)

        # Also convolve in the spatial direction.
        # FIXME: The need for this convolution might be that the simulated data is just 
        # too good, and we need to add optical aberrations.
        spat_conv_weights = np.array([0.25, .5, .25])
        if self.transpose:
            pixel_inv_var = ndimage.convolve1d(pixel_inv_var, spec_conv_weights,
                                               axis=0)
            pixel_inv_var = 1.0 / ndimage.convolve1d(1.0 / pixel_inv_var,
                                                     spat_conv_weights, axis=1)
        else:
            pixel_inv_var = ndimage.convolve1d(pixel_inv_var, spec_conv_weights,
                                               axis=1)
            pixel_inv_var = 1.0/ndimage.convolve1d(1.0/pixel_inv_var,
                                                   spat_conv_weights, axis=0)
        
        # Loop through all orders then through all y pixels.
        print("    Extracting order ", end="")
        for i in range(nm):
            print(f"{self.arm.m_min+i}...", end="")
            sys.stdout.flush()

            for j in range(ny):
                #print(f"PIXEL {self.arm.m_min+i} {j}")
                debug_this_pixel = debug_cr_pixel in [(self.arm.m_min+i, j)]

                # Check for NaNs
                if x_map[i, j] != x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue

                x_ix, phi, profiles = resample_slit_profiles_to_detector(
                    profiles, profile_y_microns, x_map[i, j] + nx // 2,
                    detpix_microns=matrices[i, j, 0, 0], debug=debug_this_pixel)
                phi /= phi.sum(axis=1)[:, np.newaxis]

                # Create our column cutout for the data and the PSF. 
                # FIXME: Is "round" most correct on the next line???
                # FIXME: Profiles should be convolved by a detector pixel,
                # but binning has to be taken into account properly!
                #x_ix = int(np.round(
                #    x_map[i, j])) - nx_cutout // 2 + \
                #       np.arange(nx_cutout, dtype=int) + nx // 2

                # CRH 20220825:  Convolve the slit to detector pixels
                # by interpolating the slit profile and integrating it 
                # across the detector pixels

                # calculate the integer number of slit profile pixels per
                # detector pixel
                #n_slit_sample = int(np.round(1.0/(profile_y_pix[1] - \
                #    profile_y_pix[0])))

                #x_ix_sample = int(np.round(
                #    x_map[i, j])) - nx_cutout // 2 + \
                #       np.arange(nx_cutout*n_slit_sample,
                #                 dtype=int)/n_slit_sample + nx // 2
                #for k in range(profiles.shape[0]):
                #    interp_prof = np.interp(
                #        x_ix_sample - x_map[i, j] - nx // 2, profile_y_pix,
                #        profiles[k])

                #    phi[k] = interp_prof.reshape(x_ix.shape[0],n_slit_sample).sum(axis=1)

                #    phi[k] /= np.sum(phi[k])
                # Deal with edge effects...
                ww = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww] = 0
                phi[:, ww] = 0.0

                # Cut out our data and inverse variance. This is where we worry
                # about whether
                # the data come with orders horizontal (default) or
                # transposed (i.e. like the
                # physical detector)
                if self.transpose:
                    col_data = data[j, x_ix]
                    col_inv_var = pixel_inv_var[j, x_ix]
                    corr = 1 if correction is None else correction[j, x_ix]
                    badpix = (np.zeros_like(col_data, dtype=bool)
                              if self.badpixmask is None else self.badpixmask[j, x_ix].astype(bool))
                else:
                    col_data = data[x_ix, j]
                    col_inv_var = pixel_inv_var[x_ix, j]                    
                    corr = 1 if correction is None else correction[x_ix, j]
                    badpix = (np.zeros_like(col_data, dtype=bool)
                              if self.badpixmask is None else self.badpixmask[x_ix, j].astype(bool))

                if debug_this_pixel:
                    profile_y_pix = profile_y_microns / matrices[i, j, 0, 0]
                    print(f"DEBUGGING {x_ix}, {j}")
                    print(f"X_MAP  {self.arm.x_map[i, j*self.arm.ybin:(j+1)*self.arm.ybin]} -> {x_map[i, j]}")
                    #print(x_map[i, j] + nx // 2 + profile_y_pix)
                    print(badpix)
                    plt.ioff()
                    fig, ax = plt.subplots()
                    ax.plot(x_ix, phi[0])
                    ax.plot(profiles[0].x / matrices[i, j, 0, 0] + x_map[i, j] + nx // 2,
                            profiles[0](profiles[0].x))
                    ax.plot(x_ix, col_data / col_data.sum())
                    plt.show()
                    plt.ion()

                # Search for additional cosmic rays here, by seeing if the data
                # look different to the model.
                col_inv_var[badpix] = 0
                additional_crs, model_amps = find_additional_crs(
                    phi, col_data, col_inv_var,
                    noise_model=noise_model, snoise=snoise, sigma=sigma,
                    reject=find_crs,
                    debug=debug_this_pixel, pixel=(self.arm.m_min+i, j))
                self.num_additional_crs += len(additional_crs)

                #Insist that variance if physical. If the model has a significant error
                #on any pixel, we don't want the weights to be completely wrong - this
                #check ensures that a reasonable extraction is possible for an imperfect
                #model. For each of the phi[object] arrays, "large" values can't
                #correspond to "small" variances. So set the minimum variance to be 
                #object shot noise.
                
                # FIXME: This is not appropriate for sky fibers, and the extraction
                # mathematics assumes that pixel variance is the same for all objects.
                # A sanity check is needed, but this isn't quite it, as 
                # minimum_variance can't be object-dependent.
                if False:
                    good_pix = np.where(col_inv_var != 0)[0]
                    
                    minimum_variance = phi.T * np.sum(phi[:, good_pix]/ \
                        col_inv_var_mat[good_pix,:], axis=0)/ \
                        np.sum(phi[:, good_pix]**2, axis=0)
                
                    col_inv_var_mat[good_pix,:] = 1./np.maximum(1./col_inv_var_mat[good_pix,:], \
                        minimum_variance[good_pix,:])

                if len(additional_crs) > 0:
                    badpix[additional_crs] = True
                    if self.badpixmask is not None:
                        if self.transpose:
                            self.badpixmask[j, x_ix[additional_crs]] |= \
                                self.cr_flag
                        else:
                            self.badpixmask[x_ix[additional_crs], j] |= \
                                self.cr_flag

                #col_inv_var_mat = np.reshape(col_inv_var.repeat(no), (x_ix.size, no))

                # Fill in the "c" matrix and "b" vector from Sharp and Birchall
                # equation 9 Simplify things by writing the sum in the
                # computation of "b" as a matrix multiplication.
                # We can do this because we're content to invert the
                # (small) matrix "c" here. 

                # For extraction weights, we only want the objects if we're
                # not sky-subtracting (even though we needed the sky for CRrej)
                #phi2 = phi[:no]
                #b_mat = phi2.T * col_inv_var_mat
                #c_mat = np.dot(phi2, phi2.T * col_inv_var_mat)

                # DEBUG - should try several columns...
                #if (i==10 and j>3000):
                #    import matplotlib.pyplot as plt
                #    plt.ion()
                #    plt.plot(col_data/np.max(col_data))
                #    plt.plot(phi[:,0]/np.max(phi[:,0]))
                #    import pdb; pdb.set_trace()
                
                # FIXME: (from Joao? Should be deleted if so)
                # Sometimes the determinant of c_mat is 0.0
                # leading to the impossibility of inverting the matrix.
                # currently the pixel weights are unmodified and those used are
                # whatever the last iteration calculation was.

                # FIXME: (from Joao? Should be deleted if so)
                # The above problem is made more complicated by the fact
                # that due to the different slit magnifications,
                # the previous weights may not have the same shape, so
                # we make copies of b_mat for pixel weights when
                # determinant is zero and shape is different.

                # If all pixels for a given object (e.g. an arc) are
                # marked as bad, c_mat can't be inverted. In this case,
                # we do the best we can with an inverse that works for
                # no cross-talk.
                #try:
                #    pixel_weights = np.dot(b_mat, np.linalg.inv(c_mat))
                #except:
                #    pixel_weights = np.dot(b_mat,
                #        np.diag(1./np.maximum(np.diag(c_mat),1e-18)))

                # FIXME: Some tilted, bright arc lines cause strange
                # weightings here... Probably OK - only strange weightings in 2D
                # really matter, and has to be re-tested once the fitted
                # spatial scale and tilt is more robust.

                pixel_weights = np.zeros((x_ix.size, no))
                phi_scaled = phi * model_amps[:, np.newaxis]
                sum_models = phi_scaled.sum(axis=0)
                col_var = noise_model(sum_models)
                col_inv_var = astrotools.divide0(1., col_var)
                col_inv_var[badpix] = 1e-18
                for ii in np.arange(no):  # for each object
                    # Avoid NaNs if there is no flux at all in a pixel
                    frac = phi_scaled[ii] / (sum_models + 1e-10)
                    if optimal:
                        pixel_weights[:, ii] = (phi[ii] * col_inv_var * frac /
                                                np.sum(phi[ii]**2 * col_inv_var))
                    else:
                        pixel_weights[:, ii] = (np.where(col_inv_var == 0, 0, frac) /
                                                np.sum(phi[ii][col_inv_var > 0]))
                    pixel_weights[badpix, ii] = 0

                # FIXME: Search here for weights that are non-zero for
                # any overlapping orders:
                if debug_this_pixel:
                    if correction is None:
                        print("SUM", np.sum(col_data))
                    else:
                        print("SUMS", np.sum(col_data), np.sum(col_data * corr))
                        print("FLATCORR VALUES")
                        print(corr)
                    print("INVERSE VARIANCE")
                    print(col_inv_var)
                    print("EXTRACTION WEIGHTS:")
                for ii, ew_one in enumerate(extraction_weights):
                    if debug_this_pixel:
                        print(ew_one[x_ix, j])
                    ew_one[x_ix, j] += pixel_weights[:, ii]
                    ew_one[x_ix, j][badpix] = 0
                    if debug_this_pixel:
                        print(ew_one[x_ix, j])

                # Actual extraction is simple: Just matrix-multiply the
                # (flat-)corrected data by the weights.
                #extracted_flux[i, j, :] = np.dot(col_data, pixel_weights)
                extracted_flux[i, j, :] = np.dot(col_data * corr, pixel_weights)

                # Rather than trying to understand and
                # document Equation 17 from Sharp and Birchall, which 
                # doesn't make a lot of sense...  lets just calculate the
                # variance in the simple explicit way for a linear combination
                # of independent pixels.
                extracted_var[i, j, :] = np.dot(
                    corr * corr / np.maximum(col_inv_var, 1e-18), pixel_weights ** 2)
                if debug_this_pixel:
                    print("EXTRACTED FLUXES", extracted_flux[i, j])
                    print("EXTRACTED VAR", extracted_var[i, j])

        print("\n")
        return extracted_flux, extracted_var, extraction_weights

    def two_d_extract(self, data=None, fl=None, extraction_weights=None,
                      vararray=None):
        """
        Perform two-dimensional flux extraction.

        Extract using 2D information. The lenslet model used is a collapsed
        profile in 1D, but we then take into account the slit shear/rotation by
        interpolating this 1D slit profile to the nearest two pixels along each
        row (y-axis in code).

        One key difference to Sharp and Birchall is that :math:`c_{kj}` (between
        equations 8 and 9) is the correct normalisation for a (fictitious)
        1-pixel wide PSF centered exactly on a pixel, but not for a continuum.
        We normalise correctly for a continuum by having one of the :math:`\phi`
        functions being one-pixel wide along the slit, and the other being
        unbounded in the dispersion direction.

        Parameters
        ----------
        data: :obj:`numpy.ndarray`, optional
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        fl: string , optional
            A fits file with conventional row/column directions containing the
            data to be extracted.
            
        extraction_weights: :obj:`numpy.ndarray`, optional
            Extraction weights created from a call to one_d_extract. Separating
            this makes the code more readable, but is not speed optimised.
        """

        # FIXME: Much of this code is doubled-up with one_d_extract.
        # Re-factoring is needed!
        if data is None:
            if fl is None:
                raise UserWarning("ERROR: Must input data or file")
            else:
                data = pyfits.getdata(fl)

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel variance if no variance is provided.
        if vararray is None:
            vararray = np.maximum(data, 0) + self.rnoise ** 2

        # Correct for scattered light - a place-holder, to show where it can
        # most easily fit.
        mask = np.sum(extraction_weights, axis=0) == 0
        data = subtract_scattered_light(data, mask)

        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            return 'Extraction failed, unable to bin models.'

        # Set up convenience local variables
        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx / self.arm.xbin)

        # Our profiles... we re-extract these in order to include the centroids.
        profile, centroids = self.slitview.slit_profile(arm=self.arm.arm, \
                                                        return_centroid=True)
        profiles = [profile]
        profile_y_microns = (np.arange(profile.size) -
                             profile.size / 2 + 0.5) * self.slitview.microns_pix

        # Convolve centroids in order to average over sub-pixel effects.
        # also convert the centroids to units of slit plane microns.
        slit_microns_per_det_pix_x = np.mean(matrices[:, :, 0, 0])
        slit_pix_per_det_pix = slit_microns_per_det_pix_x / \
                               self.slitview.microns_pix

        # Create the profile of a mean detector pixel in slit pixel units.
        det_pix = np.ones(int(slit_pix_per_det_pix) + 2)
        det_pix[0] = 0.5 * (slit_pix_per_det_pix - int(slit_pix_per_det_pix))
        det_pix[-1] = det_pix[0]
        det_pix /= np.sum(det_pix)
        for c in centroids:
            c = np.convolve(c, det_pix, mode='same') * self.slitview.microns_pix
        centroids *= self.slitview.microns_pix

        # Number of "objects"
        no = extraction_weights.shape[0]
        extracted_flux = np.zeros((nm, ny, no))
        extracted_var = np.zeros((nm, ny, no))
        extracted_covar = np.zeros((nm, ny - 1, no))

        # Mask any bad pixels with pixel inverse variance of 0.
        pixel_inv_var = 1.0 / vararray
        if self.badpixmask is not None:
            pixel_inv_var[self.badpixmask.astype(bool)] = 0.0

        # Loop through all orders then through all y pixels.
        print("    Extracting order ", end="")
        for i in range(nm):
            print(f"{self.arm.m_min+i}...", end="")
            sys.stdout.flush()

            # Create an empty weight array. Base the size on the largest slit
            # magnification for this order.
            nx_cutout = int(np.ceil(self.slitview.slit_length / np.min(
                matrices[i, :, 0, 0])))
            ny_cutout = 2 * \
                        int(nx_cutout * np.nanmax(
                            np.abs(self.slit_tilt)) / 2) + 3

            for j in range(ny):
                # Check for NaNs
                if x_map[i, j] != x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue

                # Create our cutout of columns for the data.
                x_ix = int(np.round(
                    x_map[i, j])) - nx_cutout // 2 + \
                       np.arange(nx_cutout, dtype=int) + nx // 2
                y_ix = j + np.arange(ny_cutout, dtype=int) - ny_cutout // 2

                # CJS new resampling stuff
                x_ix, _, profiles = resample_slit_profiles_to_detector(
                    profiles, profile_y_microns, x_map[i, j] + nx // 2,
                    detpix_microns = matrices[i, j, 0, 0])
                nx_cutout = x_ix.size

                # Deal with edge effects...
                ww_x = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww_x] = 0
                # phi[:, ww, :] = 0.0
                # phi1d[:, ww, :] = 0.0
                ww_y = np.where((y_ix >= ny) | (y_ix < 0))[0]
                y_ix[ww_y] = 0
                # phi[ww, :, :] = 0.0
                xy = np.meshgrid(x_ix, y_ix, indexing='ij')

                # Cut out our data and inverse variance.
                col_data = data[tuple(xy)]
                col_inv_var = pixel_inv_var[tuple(xy)]
                col_weights = np.array([ew[tuple(xy)] for ew in extraction_weights])

                col_weights[:, ww_x, :] = 0.0
                col_weights[:, :, ww_y] = 0.0
                # Find the pixel (including fractional pixels) within our
                # cutout that we'll use for extraction. First - find the pixel
                # coordinates according to slit tilt:

                ysub_pix = (x_ix - x_map[i, j] - nx // 2) * \
                            self.slit_tilt[i, j] + ny_cutout / 2 - 0.5

                # DEBUG
                # if (np.max(col_data) > 1e3):

                # Next, add the contribution of the centroid in the slit
                # viewing camera.
                # The [1,1] component of the matrix is slit_microns_per_det_
                # pix_y
                col_ix = np.arange(nx_cutout) - nx_cutout // 2

                # Interpolate onto the slit coordinates
                # FIXME: See 1d code for how this was done for profiles...
                # PRV: This is only absolutely needed for PRV mode, with 
                # matrices[i,j,1,1] coming from "specmod.fits".
                ysub_pix += np.interp(x_ix - x_map[i, j] - nx // 2,
                                      # CRH:  below was the original version but I think it's
                                      # missing a factor of the microns/pix
                                      # slit_ix / matrices[i, j, 0, 0], \
                                      profile_y_microns / matrices[i, j, 0, 0],
                                      centroids / matrices[i, j, 1, 1])

                # Make sure this is within the limits of our subarray.
                # NB We have will define ysub_pix so that the first pixel
                # covers 0 to 1, not -0.5 to +0.5 for ease of coding
                ysub_pix += 0.5
                ysub_pix = np.maximum(ysub_pix, 0)
                ysub_pix = np.minimum(ysub_pix, ny_cutout - 1 - 1e-6)

                # Create the arrays needed for interpolation.
                ysub_ix_lo = ysub_pix.astype(int)
                ysub_ix_hi = ysub_ix_lo + 1
                ysub_ix_frac = ysub_pix - ysub_ix_lo
                xsub_ix = np.arange(nx_cutout).astype(int)

                # FIXME: SERIOUS Weighting here relies on col_weights being
                # approximately correct, which it isn't for bright arc lines
                # and a tilted slit.
                # We should consider if this is "good enough" carefully.
                for k in range(no):
                    extracted_flux[i, j, k] = \
                        np.dot(col_data[xsub_ix, ysub_ix_lo],
                               col_weights[k, xsub_ix, ysub_ix_lo] * (
                                       1 - ysub_ix_frac))
                    extracted_flux[i, j, k] += \
                        np.dot(col_data[xsub_ix, ysub_ix_hi], col_weights[
                            k, xsub_ix, ysub_ix_hi] * ysub_ix_frac)
                    extracted_var[i, j, k] = np.dot(
                        (1 - ysub_ix_frac) / np.maximum(
                            col_inv_var[xsub_ix, ysub_ix_lo], 1e-18),
                        col_weights[k, xsub_ix, ysub_ix_lo] ** 2)
                    extracted_var[i, j, k] += np.dot(
                        ysub_ix_frac / np.maximum(
                            col_inv_var[xsub_ix, ysub_ix_hi], 1e-18),
                        col_weights[k, xsub_ix, ysub_ix_hi] ** 2)

        print("\n")
        return extracted_flux, extracted_var

    def find_lines(self, flux, arclines, hw=12,
                   arcfile=None, # Now dead-letter - always overridden
                   inspect=False, plots=False):
        """
        Find lines near the locations of input arc lines.
        
        This is done with Gaussian fits near the location of where lines are
        expected to be. An initial decent model must be present, likely
        the result of a manual adjustment.

        Parameters
        ----------
        flux: :obj:`numpy.ndarray`
            Flux data extracted with the 1D extractor. Just the flux, not the
            variance.

        arclines: float array
            Array containing the wavelength of the arc lines.


        arcfile: float array
            Arc file data.

        inspect: bool, optional
            If true, show display of lines and where they are predicted to fall

        plots: bool, optional
            If true, plot every gaussian fit along the way for visual inspection 

        Returns
        -------

        lines_out: float array
            Whatever used to be placed in a file.
        """
        # arcfile = flux
        # Only use the middle object.
        # In High res mode this will be the object, in std mode it's the sky
        flux = flux[:, :, 0]
        arcfile = flux
        ny = flux.shape[1]
        nm = flux.shape[0]
        nx = self.arm.szx
        lines_out = []
        # Let's try the median absolute deviation as a measure of background
        # noise if the search region is not large enough for robust median
        # determination.
        if hw < 8:
            noise_level = np.median(np.abs(flux - np.median(flux)))

        if (inspect == True) and (arcfile is None):
            raise UserWarning('Must provide an arc image for the inpection')
        if inspect or plots:
            image = np.arcsinh((arcfile - np.median(arcfile)) / 1e2)
        if inspect:
            plt.imshow(image, interpolation='nearest', aspect='auto',
                       cmap=cm.gray)

        for m_ix in range(nm):
            # Select only arc lines that should be in this order.
            filtered_arclines = arclines[
                (arclines >= self.arm.w_map[m_ix, :].min())
                & (arclines <= self.arm.w_map[m_ix, :].max())]
            # This line interpolates between filtered lines and w_map on a
            # linear array, to find the expected pixel locations of the lines.
            w_ix = np.interp(filtered_arclines, self.arm.w_map[m_ix, :],
                             np.arange(ny))
            # Ensure that lines close to the edges of the chip are not
            # considered
            ww = np.where((w_ix >= hw) & (w_ix < ny - hw))[0]
            w_ix = w_ix[ww]
            arclines_to_fit = filtered_arclines[ww]
            print('order ', m_ix)
            for i, ix in enumerate(w_ix):
                # This ensures that lines too close together are not used in the
                # fit, whilst avoiding looking at indeces that don't exist.
                if len(w_ix) > 1:
                    if (np.abs(ix - w_ix[i - 1]) < 1.5 * hw):
                        continue
                    elif i != (len(w_ix) - 1) and (
                            np.abs(ix - w_ix[i + 1]) < 1.5 * hw):
                        continue
                x = np.arange(ix - hw, ix + hw, dtype=int)
                y = flux[m_ix, x]
                # Try median absolute deviation for noise characteristics if
                # Enough pixels are available per cut.
                if hw >= 7:
                    noise_level = np.median(np.abs(y - np.median(y)))
                # Any line with peak S/N under a value is not considered.
                # e.g. 20 is considered.
                if (np.max(y) < 20 * noise_level):
                    print("Rejecting due to low SNR!")
                    continue

                g_init = models.Gaussian1D(amplitude=np.max(y), mean=x[
                    np.argmax(y)], stddev=1.5)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit_g = fitting.LevMarLSQFitter()
                    g = fit_g(g_init, x, y)
                # Wave, ypos, xpos, m, amplitude, fwhm
                xpos = nx // 2 + \
                       np.interp(g.mean.value, np.arange(ny),
                                 self.arm.x_map[m_ix])
                ypos = g.mean.value

                line_to_append = [arclines_to_fit[i], ypos, xpos, m_ix +
                                  self.arm.m_min, g.amplitude.value,
                                  g.stddev.value * 2.3548]

                # If any of the values to append are nans, don't do it and
                # just go on to the next line.
                if np.isnan(line_to_append).any():
                    continue

                # This option is here to allow the user to inspect individual
                # gaussian fits. Useful to test whether the method is working.
                if plots:
                    f, sub = plt.subplots(1, 2)
                    sub[0].plot(x, y)
                    sub[0].plot(x, g(x))
                    sub[0].axvline(ix)
                    snapshot = image[int(ix - hw * 4):int(ix + hw * 4),
                               int(xpos - 40):
                               int(xpos + 40)]
                    sub[1].imshow(np.arcsinh((snapshot - np.median(snapshot)) /
                                             1e2))
                    plt.show()

                # This part allows the user to inspect the global fit and
                # position finding. 
                if inspect:
                    plt.plot(xpos, ix, 'bx')
                    plt.plot(xpos, ypos, 'rx')
                    plt.text(xpos + 10, ypos,
                             str(arclines_to_fit[i]), color='green',
                             fontsize=10)

                lines_out.append(line_to_append)

        if inspect:
            plt.axis([0, nx, ny, 0])
            plt.show()

        return np.array(lines_out)

    def match_lines(self, all_peaks, arclines, hw=12, xcor=False, log=None):
        """
        Match peaks in data to arc lines based on proximity

        Parameters
        ----------
        all_peaks: list of :obj:`numpy.ndarray`
            Locations of peaks in each order

        arclines: float array
            Array containing the wavelength of the arc lines.

        hw: int, optional
            Number of pixels from each order end to be ignored due to proximity
            with the edge of the chip. Default was 10.

        xcor: bool
            Perform initial cross-correlation to find gross shift?

        log: logger/None

        Returns
        -------
        lines_out: 2D array
            arc line wavelength, fitted position, position along orthogonal
            direction, order, amplitude, FWHM
        """
        nx, ny = self.arm.szx, self.arm.szy
        pixels = np.arange(ny)
        lines_out = []

        for m_ix, peaks in enumerate(all_peaks):
            filtered_arclines = arclines[
                (arclines >= self.arm.w_map[m_ix].min())
                & (arclines <= self.arm.w_map[m_ix].max())]
            # This line interpolates between filtered lines and w_map on a
            # linear array, to find the expected pixel locations of the lines.
            w_ix = np.interp(filtered_arclines, self.arm.w_map[m_ix, :], pixels)
            # Ensure that lines close to the edges of the chip are not
            # considered
            ww = np.where((w_ix >= hw) & (w_ix < ny - hw))[0]
            if ww.size * len(peaks) == 0:
                continue

            w_ix = w_ix[ww]
            arclines_to_fit = filtered_arclines[ww]
            if log:
                log.stdinfo(f'order {m_ix+self.arm.m_min:2d} with '
                            f'{len(peaks):2d} peaks and {ww.size:2d} arc lines')

            # Perform cross-correlation if requested
            if xcor:
                expected = np.zeros((ny,)) + np.sum([np.exp(-(cntr - pixels) ** 2 / 18)
                                   for cntr in w_ix], axis=0)  # stddev=3
                synth = np.zeros((ny,)) + np.sum([np.exp(-(g.mean.value - pixels) ** 2 / 18)
                                for g in peaks], axis=0)  # stddev=3
                # -ve means line peaks are "ahead" of expected
                shift = np.correlate(expected, synth, mode='full')[ny-hw-1:ny+hw].argmax() - hw
            else:
                shift = 0

            xpos = lambda g: nx // 2 + np.interp(g.mean.value, pixels,
                                                 self.arm.x_map[m_ix])
            matched = matching.match_sources([g.mean.value for g in peaks],
                                             w_ix - shift, radius=hw)
            new_lines = [(arclines_to_fit[m], g.mean.value, xpos(g), m_ix + self.arm.m_min,
                          g.amplitude.value, g.stddev.value * 2.3548)
                         for i, (m, g) in enumerate(zip(matched, peaks)) if m > -1]
            lines_out.extend(new_lines)

        return np.array(lines_out)


def resample_slit_profiles_to_detector(profiles, profile_y_microns=None,
                                       profile_center=None, detpix_microns=None,
                                       debug=False):
    """
    This code takes one or more 1D profiles from the slitviewer camera and
    resamples them onto the (coarser) sampling of the echellograms. Because we
    are resampling a discretely-sampled profile onto another discrete pixel
    grid, "ringing" effects abound and so, in an attempt to reduce these, the
    profiles are expressed as cubic splines, with this function converting
    arrays into CubicSpline objects and returning those for later use (to
    prevent repeated computations; for these reason, the splines are expressed
    as functions of the offset from the centre of the slit in microns, so the
    pixel scale of the echellogram must be passed to this function).

    Parameters
    ----------
    profiles: iterable (length M) of arrays (N pixels) or spline objects
        profiles of each of the M objects, as measured on the slitviewer camera
    profile_y_microns: ndarray
        locations of each slitviewer pixel in microns, relative to slit center
    profile_center: float
        detector pixel location of center of slit (from XMOD)
    detpix_microns: float
        size of each detector pixel in microns (from SPATMOD)

    Returns
    -------
    x_ix: array of ints, of shape (P,)
        pixels in the spatial direction of the echellogram where the resampled
        profiles are to be placed
    phi: MxP ndarray
        profiles in resampled detector pixel space
    profiles: list of CubicSpline objects
        the profiles of the objects expressed as cubic splines (if they were not
        already)
    """
    # Must turn an array into a list so items can be replaced with CubicSplines
    # if they're not already
    if not isinstance(profiles, list):
        profiles = list(profiles)

    half_slitprof_pixel = 0.5 * np.diff(profile_y_microns).mean()

    # To avoid repeated computations, we recast each profile as a CubicSpline
    # with x being the separation from the slit center in *microns*. The
    # splines only extend over the region of data, to ensure they're precisely
    # zero away from this.
    for i, profile in enumerate(profiles):
        if not isinstance(profile, CubicSpline):
            # There's probably a better way to do this
            good = np.ones_like(profile, dtype=bool)
            for j, val in enumerate(profile):
                if val == 0:
                    good[j] = False
                else:
                    break
            for j, val in enumerate(profile[::-1], start=1):
                if val == 0:
                    good[-j] = False
                else:
                    break
            xspline = np.r_[profile_y_microns[good].min() - half_slitprof_pixel,
                            profile_y_microns[good],
                            profile_y_microns[good].max() + half_slitprof_pixel]
            # Apparently the antiderivative is also a CubicSpline object
            profiles[i] = CubicSpline(xspline, np.r_[0, profile[good], 0]).antiderivative()

    slitview_pixel_edges = profile_center + np.linspace(
        profile_y_microns[0] - half_slitprof_pixel,
        profile_y_microns[-1] + half_slitprof_pixel,
        profile_y_microns.size + 1) / detpix_microns
    x_ix = np.round([slitview_pixel_edges.min(),
                     slitview_pixel_edges.max() + 1]).astype(int)
    pixel_edges_microns = (np.arange(x_ix[0], x_ix[1]+0.1) - 0.5 -
                           profile_center) * detpix_microns

    phi = np.zeros((len(profiles), x_ix[1] - x_ix[0]))
    for i, profile in enumerate(profiles):
        left = np.minimum(np.maximum(pixel_edges_microns[:-1], profile.x[0]), profile.x[-1])
        boundaries = np.r_[left, max(min(pixel_edges_microns[-1], profile.x[-1]), profile.x[0])]
        phi[i] = np.diff(profile(boundaries))

    if debug:
        for j, profile in enumerate(profiles):
            print(f"PROFILE {j}")
            print(profile.x)
            print(profile(profile.x))
            print("="*60)
            #print(pixel_edges_microns)
            #print(profile(pixel_edges_microns))
            #print("="*60)
            for i in range(phi.shape[1]):
                x1 = min(max(pixel_edges_microns[i], profile.x[0]), profile.x[-1])
                x2 = max(min(pixel_edges_microns[i+1], profile.x[-1]), profile.x[0])
                print(x1, x2, profile.integrate(x1, x2))
            print(phi[j])

    phi[phi < 0] = 0
    return np.arange(*x_ix, dtype=int), phi, profiles
