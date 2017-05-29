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

class Extractor():
    def __init__(self, polyspect_instance, slitview_instance,
                 gain = 1.0, rnoise = 3.0,
                 badpixmask=[], transpose=False):
        """A class to extract data for each arm of the spectrograph.

        The extraction is defined by 3 key parameters: an "x_map", which is
        equivalent to 2dFDR's tramlines and contains a physical x-coordinate for
        every y (dispersion direction) coordinate and order, and a "w_map",
        which is the wavelength corresponding to every y (dispersion direction)
        coordinate and order.

        The slit parameters is defined by the slitview_instance. If a different
        slit profile is to be assumed, then another (e.g. simplified or fake)
        slitview instance should be made, rather than modifying this class.

        Parameters
        ----------
        polyspect_instance: Polyspect object
            This defines e.g. whether the camera is "red" or "blue"

        slitview_instance: SlitView object
            This defines whether the mode is "std" or "high"

        gain: float
            gain in electrons per ADU. From fits header.

        rnoise: float
            Expected readout noise in electrons. From fits header.

        badpixmask: numpy array
            A list of [y,x] bad pixel co-ordinates

        transpose: bool (optional)
            Do we transpose the data before extraction?
        """
        self.arm = polyspect_instance
        self.slitview = slitview_instance
        self.transpose = transpose
        self.gain = gain
        self.rnoise = rnoise
        self.badpixmask = badpixmask

        # TODO: This warning could probably be neater.
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

    def one_d_extract(self, data=[], file='', lenslet_profile='slitview',
                      correct_for_sky=True):
        """ Extract flux by integrating down columns (the "y" direction),
        using an optimal extraction method.

        Given that some of this code is in common with two_d_extract, the
        routines could easily be merged... however that would make one_d_extract
        less readable.

        Parameters
        ----------
        data: numpy array (optional)
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        file: string (optional)
            A fits file with conventional row/column directions containing the
            data to be extracted.

        correct_for_sky: bool
            Do we correct the object slit profiles for sky? Should be yes for
            objects and no for flats/arcs.

        WARNING: Binning not implemented yet"""

        if len(data) == 0:
            if len(file) == 0:
                print("ERROR: Must input data or file")
            else:
                data = pyfits.getdata(file)

        ny = self.arm.x_map.shape[1]
        nm = self.arm.x_map.shape[0]
        nx = self.arm.szx

        #Our profiles... TODO: extract x-centroids as well for PRV.
        profiles = self.slitview.object_slit_profiles(arm=self.arm.arm, correct_for_sky=correct_for_sky)

        # Number of "objects" and "slit pixels"
        no = profiles.shape[0]
        n_slitpix = profiles.shape[1]
        profile_y_microns = (np.arange(n_slitpix) -
                             n_slitpix//2)*self.slitview.microns_pix

        # To consider slit tilt for a 1D extraction, we need to know the profile
        # centroids in the "y" direction.
        y_ix = np.arange(n_slitpix) - n_slitpix//2
        y_centroids = np.empty( (no) )
        x_centroids = np.zeros( (no) )
        for y_centroid, profile in zip(y_centroids, profiles):
            y_centroid = np.sum(profile * y_ix)/np.sum(profile)
        centroids = np.array([x_centroids, y_centroids])

        #Our extracted arrays.
        extracted_flux = np.zeros((nm, ny, no))
        extracted_var = np.zeros((nm, ny, no))

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel inverse variance.
        pixel_inv_var = 1.0 / (np.maximum(data, 0) / self.gain + self.rnoise**2)
        pixel_inv_var[self.badpixmask] = 0.0

        # Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order: {0:d}".format(i))

            #Create an empty weight array. Base the size on the largest slit magnification
            #for this order.
            nx_cutout = int(np.ceil(self.slitview.slit_length/np.min(self.arm.matrices[i,:,0,0])))
            phi = np.empty((nx_cutout, no))

            for j in range(ny):
                #Compute offsets in microns directly from the y_centroid and the matrix.
                #Although this is 2D, we only worry ourselves about the
                slitim_offsets = np.dot(self.arm.matrices[i,j], centroids)
                profile_y_pix = profile_y_microns/self.arm.matrices[i,j,0,0]
                #pdb.set_trace()
                # Check for NaNs
                if self.arm.x_map[i, j] != self.arm.x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue

                # Create our column cutout for the data and the PSF. !!! Is
                # "round" correct on the next line???
                x_ix = int(np.round(
                    self.arm.x_map[i, j])) - nx_cutout // 2 +\
                    np.arange(nx_cutout, dtype=int) + nx // 2
                for k in range(no):
                    phi[:, k] = np.interp(
                        x_ix - self.arm.x_map[i, j] - nx // 2, profile_y_pix,
                        profiles[k])
                    phi[:, k] /= np.sum(phi[:, k])
                # Deal with edge effects...
                ww = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww] = 0
                phi[ww, :] = 0.0


                # Cut out our data and inverse variance. This is where we worry about whether
                #the data come with orders horizontal (default) or transposed (i.e. like the
                #physical detector)
                if self.transpose:
                    col_data = data[j, x_ix]
                    col_inv_var = pixel_inv_var[j, x_ix]
                else:
                    col_data = data[x_ix, j]
                    col_inv_var = pixel_inv_var[x_ix, j]

                # Fill in the "c" matrix and "b" vector from Sharp and Birchall
                # equation 9 Simplify things by writing the sum in the
                # computation of "b" as a matrix multiplication.
                # We can do this because we're content to invert the
                # (small) matrix "c" here. Equation 17 from Sharp and Birchall
                # doesn't make a lot of sense... so lets just calculate the
                # variance in the simple explicit way.
                col_inv_var_mat = np.reshape(
                    col_inv_var.repeat(no), (nx_cutout, no))
                b_mat = phi * col_inv_var_mat
                c_mat = np.dot(phi.T, phi * col_inv_var_mat)
                pixel_weights = np.dot(b_mat, np.linalg.inv(c_mat))
                extracted_flux[i, j, :] = np.dot(col_data, pixel_weights)
                extracted_var[i, j, :] = np.dot(
                    1.0 / np.maximum(col_inv_var, 1e-12), pixel_weights**2)

        return extracted_flux, extracted_var

    def two_d_extract(self, file='', data=[], lenslet_profile='sim', rnoise=3.0,
                      deconvolve=True):
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

        Note that the input data has to be the transpose of a conventional
        echellogram

        TODO:
        1) Neaten the approximate matrix inverse square root

        Parameters
        ----------
        data: numpy array (optional)
            Image data, transposed so that dispersion is in the "y" direction.
            Note that this is the transpose of a conventional echellogram.
            Either data or file must be given

        file: string (optional)
            A fits file with conventional row/column directions containing the
            data to be extracted.

        lenslet_profile: 'square' or 'sim'
            Shape of the profile of each fiber as used in the extraction.
            For a final implementation, 'measured' should be a possibility.
            'square' assigns each pixel uniquely to a single lenslet.
            For testing only

        rnoise: float
            The assumed readout noise.

        deconvolve: bool
            Do we deconvolve so that neighboring extracted spectral points
            are statistically independent? This is an approximate deconvolution
            (a linear function of 5 neighboring pixels) so is reasonably robust.
        """

        if len(data) == 0:
            if len(file) == 0:
                print("ERROR: Must input data or file")
            else:
                # Transpose the data from the start.
                data = pyfits.getdata(file).T

        ny = self.x_map.shape[1]
        nm = self.x_map.shape[0]
        nx = self.szx

        # Number of "objects"
        no = self.square_profile.shape[1]
        extracted_flux = np.zeros((nm, ny, no))
        extracted_var = np.zeros((nm, ny, no))
        extracted_covar = np.zeros((nm, ny - 1, no))

        # Assuming that the data are in photo-electrons, construct a simple
        # model for the pixel inverse variance.
        pixel_inv_var = 1.0 / (np.maximum(data, 0) + rnoise**2)
        pixel_inv_var[self.badpixmask] = 0.0

        # Loop through all orders then through all y pixels.
        for i in range(nm):
            print("Extracting order index: {0:d}".format(i))
            # Based on the profile we're using, create the local offsets and
            # profile vectors
            if lenslet_profile == 'sim':
                offsets = self.sim_offsets[:, i]
                profile = self.sim_profile
            else:
                print("Only sim lenslet profile available for 2D extraction so far...")
                raise userwarning
            nx_cutout = 2 * int((np.max(offsets) - np.min(offsets)) / 2) + 2
            ny_cutout = 2 * \
                int(nx_cutout * np.nanmax(np.abs(self.slit_tilt)) / 2) + 3
            for j in range(ny):
                phi = np.zeros((ny_cutout, nx_cutout, no))
                phi1d = np.zeros((ny_cutout, nx_cutout, no))
                # Check for NaNs
                if self.x_map[i, j] != self.x_map[i, j]:
                    extracted_var[i, j, :] = np.nan
                    continue
                # Create our column cutout for the data and the PSF
                x_ix = int(self.x_map[i, j]) - nx_cutout // 2 + \
                    np.arange(nx_cutout, dtype=int) + nx // 2
                y_ix = j + np.arange(ny_cutout, dtype=int) - ny_cutout // 2
                for k in range(no):
                    x_prof = np.interp(
                        x_ix - self.x_map[i, j] - nx // 2, offsets,
                        profile[:, k])
                    y_pix = (x_ix - self.x_map[i, j] - nx // 2) * \
                        self.slit_tilt[i, j] + ny_cutout // 2
                    frac_y_pix = y_pix - y_pix.astype(int)
                    subx_ix = np.arange(nx_cutout, dtype=int)
                    phi[y_pix.astype(int), subx_ix, k] = (
                        1 - frac_y_pix) * x_prof
                    phi[y_pix.astype(int) + 1, subx_ix,
                        k] = frac_y_pix * x_prof
                    phi[:, :, k] /= np.sum(phi[:, :, k])
                    x_prof /= np.sum(x_prof)
                    phi1d[:, :, k] = np.tile(x_prof, ny_cutout).reshape(
                        (ny_cutout, nx_cutout))
                # Deal with edge effects...
                ww = np.where((x_ix >= nx) | (x_ix < 0))[0]
                x_ix[ww] = 0
                phi[:, ww, :] = 0.0
                phi1d[:, ww, :] = 0.0
                ww = np.where((y_ix >= ny) | (y_ix < 0))[0]
                y_ix[ww] = 0
                phi[ww, :, :] = 0.0
                xy = np.meshgrid(y_ix, x_ix, indexing='ij')
                # Cut out our data and inverse variance.
                col_data = data[xy].flatten()
                col_inv_var = pixel_inv_var[xy].flatten()
                # Fill in the "c" matrix and "b" vector from Sharp and Birchall
                # equation 9
                # Simplify things by writing the sum in the computation of "b"
                # as a matrix
                # multiplication. We can do this because we're content to invert
                # the (small) matrix "c" here. Equation 17 from Sharp and
                # Birchall doesn't make a lot of sense...
                # so lets just calculate the variance in the simple explicit way
                col_inv_var_mat = np.reshape(
                    col_inv_var.repeat(no), (ny_cutout * nx_cutout, no))
                phi = phi.reshape((ny_cutout * nx_cutout, no))
                phi1d = phi1d.reshape((ny_cutout * nx_cutout, no))
                b_mat = phi * col_inv_var_mat
                c_mat = np.dot(phi.T, phi1d * col_inv_var_mat)
                pixel_weights = np.dot(b_mat, np.linalg.inv(c_mat))
#                if (j==1000):
#                        pdb.set_trace()
                extracted_flux[i, j, :] = np.dot(col_data, pixel_weights)
                extracted_var[i, j, :] = np.dot(
                    1.0 / np.maximum(col_inv_var, 1e-12), pixel_weights**2)
                if (j > 0):
                    extracted_covar[i, j - 1, :] = \
                        np.dot(1.0 / np.maximum(col_inv_var, 1e-12),
                               pixel_weights * np.roll(last_pixel_weights,
                                                       -nx_cutout, axis=0))
                last_pixel_weights = pixel_weights.copy()
#                if (j > 591):
#                    pdb.set_trace()
        if (deconvolve):
            # Create the diagonals of the matrix Q gradually, using the Taylor
            # approximation for the matrix inverse.
            #(Bolton and Schlegel 2009, equation 10)
            #D = diag(C)
            # A = D^{-1/2} (C-D) D^{-1/2}, so C = D^{1/2}(I + A)D^{1/2}
            # Then if Q = (I - 1/2 A + 3/8 A^2) D^{-1/2}
            #... then C^{-1} = QQ, approximately.
            # Note that all of this effort doesn't really seem to achieve much
            # at all in practice...
            # an extremely marginal improvement in resolution... but at least
            # formal pixel-to-pixel data independence is returned.
            extracted_sig = np.sqrt(extracted_var)
            a_diag_p1 = extracted_covar / \
                extracted_sig[:, :-1, :] / extracted_sig[:, 1:, :]
#            a_diag_m1 = extracted_covar/extracted_var[:,1:,:]
            Q_diag = np.ones((nm, ny, no))
            Q_diag[:, :-1, :] += 3 / 8.0 * a_diag_p1**2
            Q_diag[:, 1:, :] += 3 / 8.0 * a_diag_p1**2
#            Q_diag[:,:-1,:] += 3/8.0*a_diag_p1*a_diag_m1
#            Q_diag[:,1:,:]  += 3/8.0*a_diag_p1*a_diag_m1
            Q_diag /= extracted_sig
            extracted_sqrtsig = np.sqrt(extracted_sig)
            Q_diag_p2 = 3 / 8.0 * a_diag_p1[:, :-1, :] * a_diag_p1[
                :, 1:, :] / extracted_sqrtsig[:, 2:, :] /\
                extracted_sqrtsig[:, :-2, :]
#            Q_diag_m2 = 3/8.0*a_diag_m1[:,:-1,:]*a_diag_m1[:,1:,:]/extracted_sig[:,:-2,:]
#            Q_diag_m1 = -0.5*a_diag_m1/extracted_sig[:,:-1,:]
            Q_diag_p1 = -0.5 * a_diag_p1 / \
                extracted_sqrtsig[:, 1:, :] / extracted_sqrtsig[:, :-1, :]
    # The approximation doesn't seem to be quite right, with the ~3% uncertainty
    # on the diagonal of cinv, when there should only be a ~1% uncertainty
    # (obtained by going to the next term in the Taylor expansion).
    # But pretty close...
    # Q = np.diag(Q_diag[0,:,0]) + np.diag(Q_diag_m1[0,:,0],k=-1) +
    # np.diag(Q_diag_p1[0,:,0],k=+1) + np.diag(Q_diag_p2[0,:,0],k=+2) +
    # np.diag(Q_diag_m2[0,:,0],k=-2)
    # cinv_approx = np.dot(Q,Q)
    # cinv = np.diag(extracted_var[0,:,0]) + np.diag(extracted_covar[0,:,0],k=1)
    # + np.diag(extracted_covar[0,:,0],k=-1)
    # cinv = np.linalg.inv(cinv)
            # Now we have a sparse matrix with 5 terms. We need to sum down the
            # rows, ignoring the
            # edge pixels
#            s_vect = Q_diag[:,2:-2,:] + Q_diag_p1[:,1:-2,:] +
# Q_diag_m1[:,2:-1,:] + Q_diag_p2[:,:-2,:] + Q_diag_m2[:,2:,:]
            s_vect = Q_diag.copy()
            s_vect[:, :-1, :] += Q_diag_p1
            s_vect[:, :-2, :] += Q_diag_p2
            s_vect[:, 1:, :] += Q_diag_p1
            s_vect[:, 2:, :] += Q_diag_p2
            new_var = 1.0 / s_vect**2
            new_flux = extracted_flux * Q_diag / s_vect
            new_flux[:, :-1, :] += extracted_flux[:, 1:, :] * \
                Q_diag_p1 / s_vect[:, 1:, :]
            new_flux[:, :-2, :] += extracted_flux[:, 2:, :] * \
                Q_diag_p2 / s_vect[:, 2:, :]
            new_flux[:, 1:, :] += extracted_flux[:, :-1, :] * \
                Q_diag_p1 / s_vect[:, :-1, :]
            new_flux[:, 2:, :] += extracted_flux[:, :-2, :] * \
                Q_diag_p2 / s_vect[:, :-2, :]

            # Fill in the Variance and Flux arrays with NaNs, so that the
            # (not computed) edges are undefined.
 #           new_flux = np.empty_like(extracted_flux)
 #           new_var = np.empty_like(extracted_var)
 #           new_flux[:,:,:]=np.nan
 #           new_var[:,:,:]=np.nan
            # Now fill in the arrays.
 #           new_var[:,2:-2,:] = 1.0/s_vect**2
 #           new_flux[:,2:-2,:] =  extracted_flux[:,2:-2,:]*Q_diag[:,2:-2,:]/s_vect
            #
 #           new_flux[:,2:-2,:] += extracted_flux[:,1:-3,:]*Q_diag_p1[:,1:-2,:]/s_vect
 #           new_flux[:,2:-2,:] += extracted_flux[:,3:-1,:]*Q_diag_p1[:,2:-1,:]/s_vect
 #           new_flux[:,2:-2,:] += extracted_flux[:,:-4,:] *Q_diag_p2[:,:-2,:]/s_vect
 #           new_flux[:,2:-2,:] += extracted_flux[:,4:,:]  *Q_diag_p2[:,2:,:]/s_vect

            return new_flux, new_var
        else:
            return extracted_flux, extracted_var

    def find_lines(self, flux, arclines, hw=10, arcfile=[],
                   inspect=False):
        """Find lines near the locations of input arc lines.

        Parameters
        ----------
        flux: numpy array
            Flux data extracted with the 1D extractor. Just the flux, not the
            variance.

        arclines: float array
            Array containing the wavelength of the arc lines.

        hw: int (optional)
            Number of pixels from each order end to be ignored due to proximity
            with the edge of the chip.

        flat_data: float array
            Flat field data? For what?

        arcfile: float array
            Arc file data.

        inspect: boolean (optional)
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

        if inspect and len(arcfile) == 0:
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
                #pdb.set_trace()
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
