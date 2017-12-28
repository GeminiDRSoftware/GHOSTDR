"""Given (x,wave,matrices, slit_profile), extract the flux from each order.
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pdb
from astropy.modeling import models, fitting
import matplotlib.cm as cm
import warnings
import scipy.ndimage as ndimage

def find_additional_crs(phi, slitim_offsets, col_data, col_inv_var, \
    snoise=0.1, nsigma=6, debug=False):
    """Utility function to search for additional cosmic rays.
    
    Parameters
    ----------
    phi: :obj:`numpy.ndarray`
        A (npix x nobj) model PSF array.
    col_data: :obj:`numpy.ndarray'
        A (npix) data array
    col_inv_var: :obj:`numpy.ndarray'
        A (npix) data variance array
    multiplicative_std: float
        An additional error variance, due to model uncertainties, which 
        are added to the data
        
    Returns
    -------
    new_bad: :obj:`numpy.ndarray'
        List of bad pixels.
    """
    n_o = phi.shape[1] #Number of objects.
    n_x = len(col_data)
    
    var_use = 1/col_inv_var + (snoise * col_data)**2

    #Create a model matrix for linear regression.
    obj_centers = slitim_offsets[1] + n_x//2
    x_ix = np.arange(n_x)
    x_mat = np.empty( (2*n_o, n_x) )
    for o_ix in range(n_o):
        x_mat[o_ix] = phi[:, o_ix]
        x_mat[o_ix + n_o] = phi[:, o_ix] * (x_ix - obj_centers[o_ix])
        
    # Now we fit a model to the col_data with var_use, using standard
    # linear regression. FIXME: Surely there is a reasonably scipy helper
    # function???
    x_mat = x_mat.T
    #FIXME: Do we use weights? If so, faster computation is possible as W is sparse.
    # performance is improved with no weights.
    # w_mat = np.diag(1./var_use)
    # xtw = np.dot(x_mat.T, w_mat) 
    #beta = np.dot(np.dot(np.linalg.inv(np.dot(xtw, x_mat)), xtw), col_data)
    beta = np.dot(np.dot(np.linalg.inv(np.dot(x_mat.T, x_mat)), x_mat.T), col_data)
    y_hat = np.maximum(np.dot(x_mat, beta),0)
    
    new_bad = np.where(col_data > y_hat + nsigma * np.sqrt(var_use))[0]
    if debug and len(new_bad)>0:
        plt.clf()
        plt.plot(y_hat, label='exp')
        plt.plot(col_data, label='data')
        plt.plot(y_hat + nsigma*np.sqrt(var_use), label='limit')
        plt.pause(.001)
    
    return new_bad

def subtract_scattered_light(data, mask):
    """Subtract a simple linear approximation to the scattered light
    on the detector. This is a place-holder function, which can be upgraded to 
    higher order scattered light subtraction, depending on the performance of 
    the spectrograph.
    
    Parameters
    ----------
    data: :obj:`numpy.ndarray'
        The data to be corrected for scattered light.
    
    mask: :obj:`numpy.ndarray'
        A mask which is True wherever there is only scattered light and no data.
    
    Returns
    -------
    data: :obj:`numpy.ndarray'
        A scattered-light corrected image.
    
    """
    return data


class Extractor():
    def __init__(self, polyspect_instance, slitview_instance,
                 gain = 1.0, rnoise = 3.0, cr_flag = 8,
                 badpixmask=None, transpose=False,
                 vararray=None):
        """A class to extract data for each arm of the spectrograph.

        The extraction is defined by 3 key parameters: an "x_map", which is
        equivalent to 2dFDR's tramlines and contains a physical x-coordinate for
        every y (dispersion direction) coordinate and order, and a "w_map",
        which is the wavelength corresponding to every y (dispersion direction)
        coordinate and order.

        The slit parameters is defined by the slitview_instance. If a different
        slit profile is to be assumed, then another (e.g. simplified or fake)
        slitview instance should be made, rather than modifying this class.

        Attributes
        ----------
        arm: Polyspect object
            This defines e.g. whether the camera is "red" or "blue"

        slitview: SlitView object
            This defines whether the mode is "std" or "high"

        gain: float, optional
            gain in electrons per ADU. From fits header.

        rnoise: float, optional
            Expected readout noise in electrons. From fits header.

        badpixmask: :obj:`numpy.ndarray', optional
            A data quality plane, which evaluates to False (i.e. 0) for good
            pixels.

        vararray: :obj:`numpy.ndarray', optional
            A variance array. 
        transpose: bool , optional
            Do we transpose the data before extraction?
        
        cr_flag: integer, optional
            When we flag additional cosmic rays in the badpixmask, what value
            should we use?
        """
        self.arm = polyspect_instance
        self.slitview = slitview_instance
        self.transpose = transpose
        self.gain = gain
        self.rnoise = rnoise
        self.vararray = vararray
        self.badpixmask = badpixmask
        self.cr_flag = cr_flag
        
        # FIXME: This warning could probably be neater.
        if not isinstance(self.arm.x_map,np.ndarray):
            raise UserWarning('Input polyspect_instance requires \
            spectral_format_with matrix to be run.')

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
        """ Function used to artificially bin the models so that they apply to
        whatever binning mode the data are. This requires knowledge of the 
        x and y binning from the arm class, which is assumed this class
        inherits. 
        
        The binning is done as essentially a running average, in which the 
        values for each binned pixel are assumed to be equivalent to the average
        value of all physical pixels that are part of the binned pixel.

        Parameters
        ----------
            None, assumes all necessary parameters are class attributes. 
        Returns
        -------
        x_map: :obj:`numpy.ndarray'
            Binned version of the x map model
        wmap: :obj:`numpy.ndarray'
            Binned version of the wavelength map model
        blaze: :obj:`numpy.ndarray'
            Binned version of the blaze map model
        matrices: :obj:`numpy.ndarray'
            Binned version of the matrices array
        """
        if self.arm.xbin==1 and self.arm.ybin==1:
            return self.arm.x_map, self.arm.w_map, self.arm.blaze,\
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
                                          int(self.arm.szy/self.arm.ybin),
                                          self.arm.ybin),axis=2)

            # Now do the same for the wavelength scale and blaze where necessary
            w_map = np.mean(w_map.reshape(n_orders,
                                          int(self.arm.szy/self.arm.ybin),
                                          self.arm.ybin),axis=2)

            blaze = np.mean(blaze.reshape(n_orders,
                                          int(self.arm.szy/self.arm.ybin),
                                          self.arm.ybin),axis=2)
            # The matrices are a bit harder to work with, but still the same
            # principle applies.
            matrices = np.mean(matrices.reshape(n_orders,
                                                int(self.arm.szy/self.arm.ybin),
                                                self.arm.ybin, 2, 2),axis=2)

        if self.arm.xbin > 1:
            # Now, naturally, the actualy x values must change according to the
            # xbin
            x_map /= self.arm.xbin
            # The w_map and blaze remain unchanged by this. 

        # Now we must modify the values of the [0,0] and [1,1] elements of
        # each matrix according to the binning to reflect the now size of
        # binned pixels.
        rescale_mat = np.array([[self.arm.xbin, 0] , [0 , self.arm.ybin]])
        matrices = np.dot(matrices, rescale_mat)

        return x_map, w_map, blaze, matrices    

            
    def one_d_extract(self, data=None, file=None, correct_for_sky=True,\
        debug_crs=False):
        """ Extract flux by integrating down columns (the "y" direction),
        using an optimal extraction method.

        Given that some of this code is in common with two_d_extract, the
        routines could easily be merged... however that would make one_d_extract
        less readable.

        Parameters
        ----------
        data: :obj:`numpy.ndarray' , optional
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        file: string , optional
            A fits file with conventional row/column directions containing the
            data to be extracted.

        correct_for_sky: bool, optional
            Do we correct the object slit profiles for sky? Should be yes for
            objects and no for flats/arcs.

        WARNING: Binning not implemented yet

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

        # TODO: implement binning
        
        if data is None:
            if file is None:
                raise UserWarning("Must input data or file")
            else:
                data = pyfits.getdata(file)

        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            return 'Extraction failed, unable to bin models.'

        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx/self.arm.xbin)

        #Our profiles... 
        #FIXME: Consider carefully whether there is a way to extract x-centroids
        #as well for PRV, as part of slitim_offsets below.
        profiles = self.slitview.object_slit_profiles(\
            arm = self.arm.arm, correct_for_sky = correct_for_sky)

        # Number of "objects" and "slit pixels"
        no = profiles.shape[0]
        n_slitpix = profiles.shape[1]
        profile_y_microns = (np.arange(n_slitpix) -
                             n_slitpix//2)*self.slitview.microns_pix

        
        # To consider slit tilt for a 1D extraction, we need to know the profile
        # centroids in the "y" direction. i.e. In principle, we could modify the
        # wavelength scale for each object based on this. If 2D extraction works
        # well, such an approach is not needed, but lets keep the idea of this
        # code here for now.
        
        # FIXME: This part doesn't actually do anything. But it's also not used.
        
        #y_ix = np.arange(n_slitpix) - n_slitpix//2
        y_centroids = np.empty( (no) )
        x_centroids = np.zeros( (no) )
        for object_ix, profile in enumerate(profiles):
            y_centroids[object_ix] = np.sum(profile *
                                            profile_y_microns)/np.sum(profile)
        centroids = np.array([x_centroids, y_centroids])

        #Our extracted arrays, and the weights array
        extracted_flux = np.zeros((nm, ny, no))
        extracted_var = np.zeros((nm, ny, no))
        extraction_weights = np.zeros((no, nx, ny))

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel inverse variance.
        # This really should come from an input "vararray" because of differing
        # gain and readout noise parameters for different amplifiers, and known
        # bad pixels.
        if self.vararray is None:
            pixel_inv_var = 1.0 / (np.maximum(data, 0) /
                                   self.gain + self.rnoise**2)
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
        # Also convolve in the spatial direction, to 
        # FIXME: The need for this convolution might be that the simulated data is just 
        # too good, and we need to add optical aberrations.
        spat_conv_weights = np.array([0.25,.5,.25])
        if self.transpose:
            pixel_inv_var = ndimage.convolve1d(pixel_inv_var, spec_conv_weights,
                                               axis=0)
            pixel_inv_var = 1.0/ndimage.convolve1d(1.0/pixel_inv_var,
                                                   spat_conv_weights, axis=1)
        else:
            pixel_inv_var = ndimage.convolve1d(pixel_inv_var, spec_conv_weights,
                                               axis=1)
            pixel_inv_var = 1.0/ndimage.convolve1d(1.0/pixel_inv_var,
                                                   spat_conv_weights, axis=0)
        
        if self.badpixmask is not None:
            pixel_inv_var[self.badpixmask.astype(bool)] = 0.1
            
        # Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order: {0:d}".format(i))

            #Create an empty weight array. Base the size on the largest slit magnification
            #for this order.
            nx_cutout = int(np.ceil(self.slitview.slit_length/np.min(matrices[i,:,0,0])))
            phi = np.empty((nx_cutout, no))

            for j in range(ny):
                # Compute offsets in slit-plane microns directly from the y_centroid and 
                # the matrix.
                # FIXME: It is unclear what to DO with this for a 1D extraction, unless
                # we were going to output a new wavelength scale associated with a
                # 1D extraction. This is currently only used to make the fitting as 
                # part of the CR rejection neat.
                slitim_offsets = np.dot(np.linalg.inv(matrices[i,j]),
                                        centroids)
                
                profile_y_pix = profile_y_microns/matrices[i,j,0,0]
                
                # Check for NaNs
                if x_map[i, j] != x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue

                # Create our column cutout for the data and the PSF. 
                # FIXME: Is "round" most correct on the next line???
                # FIXME: Profiles should be convolved by a detector pixel,
                # but binning has to be taken into account properly!
                x_ix = int(np.round(
                    x_map[i, j])) - nx_cutout // 2 +\
                    np.arange(nx_cutout, dtype=int) + nx // 2
                for k in range(no):
                    phi[:, k] = np.interp(
                        x_ix - x_map[i, j] - nx // 2, profile_y_pix,
                        profiles[k])
                    phi[:, k] /= np.sum(phi[:, k])
                # Deal with edge effects...
                ww = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww] = 0
                phi[ww, :] = 0.0

                # Cut out our data and inverse variance. This is where we worry about whether
                # the data come with orders horizontal (default) or transposed (i.e. like the
                # physical detector)
                if self.transpose:
                    col_data = data[j, x_ix]
                    col_inv_var = pixel_inv_var[j, x_ix]
                else:
                    col_data = data[x_ix, j]
                    col_inv_var = pixel_inv_var[x_ix, j]

                # Search for additional cosmic rays here, by seeing if the data
                # look different to the model.
                additional_crs = find_additional_crs(phi, slitim_offsets, col_data, \
                    col_inv_var, debug=debug_crs)
                if (len(additional_crs) > 0):
                    col_inv_var[additional_crs] = 0
                    if self.badpixmask is not None:
                        if self.transpose:
                            self.badpixmask[j, x_ix[additional_crs]] |= self.cr_flag
                        else:
                            self.badpixmask[x_ix[additional_crs], j] |= self.cr_flag

                # Fill in the "c" matrix and "b" vector from Sharp and Birchall
                # equation 9 Simplify things by writing the sum in the
                # computation of "b" as a matrix multiplication.
                # We can do this because we're content to invert the
                # (small) matrix "c" here. R
                # FIXME: Transpose matrices appropriately so that multiplication
                # order is the same as documentation.
                col_inv_var_mat = np.reshape(
                    col_inv_var.repeat(no), (nx_cutout, no))
                b_mat = phi * col_inv_var_mat
                c_mat = np.dot(phi.T, phi * col_inv_var_mat)
                # FIXME: Sometimes the determinant of c_mat is 0.0
                # leading to the impossibility of inverting the matrix.
                # currently the pixel weights are unmodified and those used are
                # whatever the last iteration calculation was.
                # Naturally, if this happens on iteration one, this will crash
                # on the next loop.
                try: pixel_weights = np.dot(b_mat, np.linalg.inv(c_mat))
                except:
                    try:
                        pixel_weights
                    except:
                        pixel_weights = b_mat.copy()
                    else:
                        pass

                #FIXME Bugshooting: Some tilted, bright arc lines cause strange
                #weightings here... Probably OK - only strange weightings in 2D
                #really matter, and has to be re-tested once the fitted 
                #spatial scale and tilt is more robust.
                
                #DEBUG
                #if (np.max(pixel_weights) > 2) and (j > 0.1*ny):
                #if (j==2000):
                #    import pdb; pdb.set_trace()
                
                #FIXME: Search here for weights that are non-zero for overlapping
                #orders
                for ii, ew_one  in enumerate(extraction_weights):
                    ew_one[x_ix, j] += pixel_weights[:,ii]

                #Actual extraction is simple: Just matrix-multiply the data by the
                #weights. 
                extracted_flux[i, j, :] = np.dot(col_data, pixel_weights)

                # Rather than trying to understand and
                # document Equation 17 from Sharp and Birchall, which 
                # doesn't make a lot of sense... so lets just calculate the
                # variance in the simple explicit way for a linear combination of
                # independent pixels.
                extracted_var[i, j, :] = np.dot(
                    1.0 / np.maximum(col_inv_var, 1e-12), pixel_weights**2)

        return extracted_flux, extracted_var, extraction_weights

    def two_d_extract(self, data=None, file=None, extraction_weights=None):
        """ Extract using 2D information. The lenslet model used is a collapsed
        profile, in 1D but where we take into account the slit shear/rotation by
        interpolating this 1D slit profile to the nearest two pixels along each
        row (y-axis in code).

        One key difference to Sharp and Birchall is that c_kj (between equations
        8 and 9) is the correct normalisation for a (fictitious) 1-pixel wide
        PSF centered exactly on a pixel, but not for a continuum. We normalise
        correctly for a continuum by having one of the \phi functions being
        one-pixel wide along the slit, and the other being unbounded in the
        dispersion direction.

        Parameters
        ----------
        data: :obj:`numpy.ndarray', optional
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        file: string , optional
            A fits file with conventional row/column directions containing the
            data to be extracted.
            
        extraction_weights: :obj:`numpy.ndarray', optional
            Extraction weights created from a call to one_d_extract. Separating
            this makes the code more readable, but is not speed optimised.
        """
            
        #FIXME: Much of this code is doubled-up with one_d_extract. Re-factoring
        #is needed!
        if data is None:
            if file is None:
                raise UserWarning("ERROR: Must input data or file")
            else:
                data = pyfits.getdata(file)
        
        # Correct for scattered light - a place-holder, to show where it can 
        # most easily fit.
        mask = np.sum(extraction_weights,axis=0) == 0
        data = subtract_scattered_light(data, mask)

        try:
            x_map, w_map, blaze, matrices = self.bin_models()
        except Exception:
            return 'Extraction failed, unable to bin models.'

        #Set up convenience local variables
        ny = x_map.shape[1]
        nm = x_map.shape[0]
        nx = int(self.arm.szx/self.arm.xbin)

        # Our profiles... we re-extract these in order to include the centroids.
        profile, centroids = self.slitview.slit_profile(arm=self.arm.arm, \
            return_centroid=True)

        # Convolve centroids in order to average over sub-pixel effects.
        # also convert the centroids to units of slit plane microns.
        slit_microns_per_det_pix_x = np.mean(matrices[:, :, 0, 0])
        slit_pix_per_det_pix = slit_microns_per_det_pix_x/self.slitview.microns_pix
        #Create the profile of a mean detector pixel in slit pixel units.
        det_pix = np.ones(int(slit_pix_per_det_pix) + 2)
        det_pix[0] = 0.5*(slit_pix_per_det_pix - int(slit_pix_per_det_pix))
        det_pix[-1] = det_pix[0]
        det_pix /= np.sum(det_pix)
        for c in centroids:
            c = np.convolve(c, det_pix, mode='same')*self.slitview.microns_pix

        # Number of "objects"
        no = extraction_weights.shape[0]
        extracted_flux = np.zeros((nm, ny, no))
        extracted_var = np.zeros((nm, ny, no))
        extracted_covar = np.zeros((nm, ny - 1, no))

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel inverse variance.
        pixel_inv_var = 1.0 / (np.maximum(data, 0) + self.rnoise**2)
        if self.badpixmask is not None:
            pixel_inv_var[self.badpixmask.astype(bool)] = 0.0
        
        # Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order: {0:d}".format(i))

            #Create an empty weight array. Base the size on the largest slit magnification
            #for this order.
            nx_cutout = int(np.ceil(self.slitview.slit_length/np.min(matrices[i,:,0,0])))
            ny_cutout = 2 * \
                int(nx_cutout * np.nanmax(np.abs(self.slit_tilt)) / 2) + 3

            for j in range(ny):
                # Check for NaNs
                if x_map[i, j] != x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue

                # Create our cutout of columns for the data. 
                x_ix = int(np.round(
                    x_map[i, j])) - nx_cutout // 2 +\
                    np.arange(nx_cutout, dtype=int) + nx // 2
                y_ix = j + np.arange(ny_cutout, dtype=int) - ny_cutout // 2

                # Deal with edge effects...
                ww = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww] = 0
                #phi[:, ww, :] = 0.0
                #phi1d[:, ww, :] = 0.0
                ww = np.where((y_ix >= ny) | (y_ix < 0))[0]
                y_ix[ww] = 0
                #phi[ww, :, :] = 0.0
                xy = np.meshgrid(x_ix, y_ix, indexing='ij')
                
                # Cut out our data and inverse variance.
                col_data = data[xy]
                col_inv_var = pixel_inv_var[xy]
                col_weights = np.array([ew[xy] for ew in extraction_weights])
                
                # Find the pixel (including fractional pixels) within our cutout that 
                # we'll use for extraction. First - find the pixel coordinates 
                # according to slit tilt:
                
                ysub_pix = (x_ix - x_map[i, j] - nx // 2) * \
                        self.slit_tilt[i, j] + ny_cutout // 2
                        
                #DEBUG
                #if (np.max(col_data) > 1e3): 
                #    import pdb; pdb.set_trace()

                # Next, add the contribution of the centroid in the slit viewing camera.
                # The [1,1] component of the matrix is slit_microns_per_det_pix_y 
                # Above, slit_ix was called profile_y_pix.
                slit_ix = np.arange(len(centroids)) - len(centroids)//2
                col_ix = np.arange(nx_cutout) - nx_cutout//2
                
                # Interpolate onto the slit coordinates
                # FIXME: See 1d code for how this was done for profiles...
                # PRV: This is only absolutely needed for PRV mode, with 
                # matrices[i,j,1,1] coming from "specmod.fits".
                ysub_pix += np.interp(x_ix - x_map[i, j] - nx // 2, \
                    slit_ix/matrices[i, j, 0, 0], \
                    centroids/matrices[i, j, 1, 1])
                
                #Make sure this is within the limits of our subarray.
                ysub_pix = np.maximum(ysub_pix, 0)
                ysub_pix = np.minimum(ysub_pix, ny_cutout - 1e-6)

                #Create the arrays needed for interpolation.
                ysub_ix_lo = ysub_pix.astype(int)
                ysub_ix_hi = ysub_ix_lo + 1
                ysub_ix_frac = ysub_pix - ysub_ix_lo 
                xsub_ix = np.arange(nx_cutout).astype(int)

                # FIXME: SERIOUS Weighting here relies on col_weights being approximately 
                # correct, which it isn't for bright arc lines and a tilted slit.
                # We should consider if this is "good enough" carefully.
                for k in range(no):
                    extracted_flux[i, j, k] = \
                        np.dot(col_data[xsub_ix, ysub_ix_lo], col_weights[k, xsub_ix, ysub_ix_lo]*(1 - ysub_ix_frac))
                    extracted_flux[i, j, k] += \
                        np.dot(col_data[xsub_ix, ysub_ix_hi], col_weights[k, xsub_ix, ysub_ix_hi] * ysub_ix_frac)
                    extracted_var[i, j, k] = np.dot(
                        (1 - ysub_ix_frac) / np.maximum(col_inv_var[xsub_ix, ysub_ix_lo], 1e-12), \
                        col_weights[k, xsub_ix, ysub_ix_lo]**2)
                    extracted_var[i, j, k] += np.dot(
                        ysub_ix_frac / np.maximum(col_inv_var[xsub_ix, ysub_ix_hi], 1e-12), \
                        col_weights[k, xsub_ix, ysub_ix_hi]**2)
                
        return extracted_flux, extracted_var     

    def find_lines(self, flux, arclines, hw=10, arcfile=None,
                   inspect=False):
        """Find lines near the locations of input arc lines. 
        
        This is done with Gaussian fits near the location of where lines are
        expected to be. An initial decent model must be present, likely
        the result of a manual adjustment.

        Parameters
        ----------
        flux: :obj:`numpy.ndarray'
            Flux data extracted with the 1D extractor. Just the flux, not the
            variance.

        arclines: float array
            Array containing the wavelength of the arc lines.

        hw: int, optional
            Number of pixels from each order end to be ignored due to proximity
            with the edge of the chip.

        arcfile: float array
            Arc file data.

        inspect: boolean, optional
            If true, show display of lines and where they are predicted to fall

        Returns
        -------

        lines_out: float array
            Whatever used to be placed in a file.
        """

        # Only use the middle object.
        # In High res mode this will be the object, in std mode it's the sky
        flux = flux[:, :, 2]
        ny = flux.shape[1]
        nm = flux.shape[0]
        nx = self.arm.szx
        lines_out = []
        # Let's try the median absolute deviation as a measure of background
        # noise.
        noise_level = np.median(np.abs(flux-np.median(flux)))

        if (inspect == True) and (arcfile is None):
            print('Must provide an arc image for the inpection')
            raise UserWarning
        if inspect:
            plt.imshow(np.arcsinh((arcfile - np.median(arcfile)) /
                                  1e2),
                       interpolation='nearest', aspect='auto', cmap=cm.gray)


        for m_ix in range(nm):
            filtered_arclines=arclines[(arclines >=  self.arm.w_map[m_ix, :].min())
                                       & (arclines <= self.arm.w_map[m_ix, :].max())]
            w_ix = np.interp(filtered_arclines, self.arm.w_map[m_ix, :], np.arange(ny))
            ww = np.where((w_ix >= hw) & (w_ix < ny - hw))[0]
            w_ix = w_ix[ww]
            arclines_to_fit = filtered_arclines[ww]
            print('order ', m_ix)
            for i, ix in enumerate(w_ix):
                # This ensures that lines too close together are not used in the
                # fit, whilst avoiding looking at indeces that don't exist.
                if (np.abs(ix-w_ix[i-1])<1.3*hw):
                    continue
                elif i!=(len(w_ix)-1) and (np.abs(ix-w_ix[i+1])<1.3*hw):
                    continue
                x = np.arange(ix - hw, ix + hw, dtype=np.int)
                y = flux[m_ix, x]
                # Any line with peak S/N under a value is not considered.
                # And reject any saturated lines.
                if (np.max(y) < 6 * noise_level) or (np.max(y) > 6E4):
                    continue
                g_init = models.Gaussian1D(amplitude=np.max(y), mean=x[
                                           np.argmax(y)], stddev=1.5)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit_g = fitting.LevMarLSQFitter()
                    g = fit_g(g_init, x, y)
                #Wave, ypos, xpos, m, amplitude, fwhm
                xpos = nx // 2 + \
                    np.interp(g.mean.value, np.arange(ny), self.arm.x_map[m_ix])
                ypos = g.mean.value
                if inspect:
                    plt.plot(xpos, ix, 'bx')
                    plt.plot(xpos, ypos, 'rx')
                    plt.text(xpos + 10, ypos,
                             str(arclines_to_fit[i]), color='green', fontsize=10)
                lines_out.append([arclines_to_fit[i], ypos, xpos, m_ix +
                                  self.arm.m_min, g.amplitude.value,
                                  g.stddev.value * 2.3548])
        if inspect:
            plt.axis([0, nx, ny, 0])
            plt.show()
        return np.array(lines_out)

