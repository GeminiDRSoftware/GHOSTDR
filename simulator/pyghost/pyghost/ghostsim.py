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

try:
    import pyfits as pf
except ImportError:
    import astropy.io.fits as pf

import pyghost.optics as optics
import pyghost.cosmic as cosmic

LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))


class Arm(object):
    """A class for each arm of the spectrograph. The initialisation function
    takes a single string representing the configuration. For GHOST, it can be
    "red" or "blue"."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, arm):
        self.arm = arm
        self.d_y = 1000/52.67          # Distance in microns
        self.theta = 65.0              # Blaze angle
        self.assym = 1.0/0.41  # Magnification
        self.gamma = 0.56      # Echelle gamma
        self.nwave = 1e2       # Wavelengths per order for interpolation.
        self.f_col = 1750.6    # Collimator focal length.
        self.lenslet_high_size = 118.0  # Lenslet flat-to-flat in microns
        self.lenslet_std_size = 197.0   # Lenslet flat-to-flat in microns
        self.microns_pix = 2.0  # slit image microns per pixel
        self.microns_arcsec = 400.0  # slit image plane microns per arcsec
        self.im_slit_sz = 2048  # Size of the image slit size in pixels.
        if arm == 'red':
            # Additional slit rotation across an order needed to match Zemax.
            self.extra_rot = 3.0
            self.szx = 6144
            self.szy = 6160
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.drot = -2.0       # Detector rotation
            self.d_x = 1000/565.   # VPH line spacing
            self.theta_i = 30.0    # Prism incidence angle
            self.alpha1 = 0.0      # First prism apex angle
            self.alpha2 = 0.0      # Second prism apex angle
            self.order_min = 34
            self.order_max = 67
        elif arm == 'blue':
            # Additional slit rotation accross an order needed to match Zemax.
            self.extra_rot = 2.0
            self.szx = 4096
            self.szy = 4112
            self.f_cam = 264.0
            self.px_sz = 15e-3
            self.d_x = 1000/1137.  # VPH line spacing
            self.theta_i = 30.0    # Prism incidence angle
            self.drot = -2.0       # Detector rotation.
            self.alpha1 = 0.0      # First prism apex angle
            self.alpha2 = 0.0      # Second prism apex angle
            self.order_min = 63
            self.order_max = 95
        else:
            print("Unknown spectrograph arm!")
            raise UserWarning

    def spectral_format(self, xoff=0.0, yoff=0.0, ccd_centre=None):
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
        orders = np.arange(self.order_min, self.order_max+1)
        wave_mins = 2*self.d_y*np.sin(np.radians(self.theta))/(orders + 1.0)
        wave_maxs = 2*self.d_y*np.sin(np.radians(self.theta))/(orders - 1.0)
        wave = np.empty((len(orders), self.nwave))
        for i in range(len(orders)):
            wave[i, :] = np.linspace(wave_mins[i], wave_maxs[i], self.nwave)
        wave = wave.flatten()
        orders = np.repeat(orders, self.nwave)
        order_frac = \
            np.abs(orders - 2*self.d_y*np.sin(np.radians(self.theta))/wave)
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
            w_ix = np.where((ypx < self.szy/2) * (ypx > -self.szy/2))[0]
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
            ypx_min = np.max([np.min(ypx[w_ix]).astype(int), -self.szy/2])
            ypx_max = np.min([np.max(ypx[w_ix]).astype(int), self.szy/2])
            y_int_m = np.arange(ypx_min, ypx_max, dtype=int)
            y_ix = y_int_m + self.szy/2
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
            #w_ix = np.where(isbad[i, :] == False)[0]
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
                    matrices[i, j, :, :] = np.zeros_like(amat)
                else:
                    matrices[i, j, :, :] = np.linalg.inv(amat)
        return x_c, w_c, b_c, matrices

    def make_lenslets(self, fluxes=[], mode='', seeing=0.8, llet_offset=0):
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
        s32 = np.sqrt(3)/2
        hex_scale = 1.15
        # equivalent to a 1 degree FWHM for an f/3 input ???
        # !!! Double-check !!!
        conv_fwhm = 30.0
        if len(fluxes) == 28:
            mode = 'high'
        elif len(fluxes) == 17:
            mode = 'std'
        elif len(mode) == 0:
            print("Error: 17 or 28 lenslets needed... or mode should be set")
            raise UserWarning
        if mode == 'std':
            n_lenslets = 17
            lenslet_width = self.lenslet_std_size
            yoffset = \
                (lenslet_width/self.microns_pix/hex_scale *
                 np.array([0, -s32, s32, 0, -s32, s32, 0])).astype(int)
            xoffset = \
                (lenslet_width/self.microns_pix/hex_scale *
                 np.array([-1, -0.5, -0.5, 0, 0.5, 0.5, 1.0])).astype(int)
        elif mode == 'high':
            n_lenslets = 28
            lenslet_width = self.lenslet_high_size
            #Order is 6 from the outer array then inner 7 then 6 from the
            #outer array. Should confirm to AAO.
            yoffset = (lenslet_width/self.microns_pix/hex_scale*s32 *
                       np.array([-2, 2, -2, -1, -1, 0,
                                 -1, -1, 0, 0, 0, 1, 1,
                                 0, 1, 1, 2, -2, 2])).astype(int)
            xoffset = (lenslet_width/self.microns_pix/hex_scale*0.5 *
                       np.array([-2, 0, 2, -3, 3, -4,
                                 -1, 1, -2, 0, 2, -1, 1,
                                 4, -3, 3, -2, 0, 2])).astype(int)
        else:
            print("Error: mode must be standard or high")

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
        h_long[:, szx/2-szy/2:szx/2+szy/2] = h_cutout
        h_long_big[:, szx/2-szy/2:szx/2+szy/2] = hbig_cutout
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
            im_cutout = optics.moffat2d(szy, \
                seeing * self.microns_arcsec / self.microns_pix / 2, beta=4.0)
            im_object[:, szx/2-szy/2:szx/2+szy/2] = im_cutout
            fluxes = np.ones(len(xoffset))

        # Go through the flux vector and fill in each lenslet.
        for i, flux in enumerate(fluxes):
            im_one = np.zeros((szy, szx))
            im_cutout = np.roll(np.roll(im_object, yoffset[i], axis=0),
                                xoffset[i], axis=1) * h_long
            im_cutout = im_cutout[szy/2-cutout_hw:szy/2+cutout_hw,
                                  szx/2-cutout_hw:szx/2+cutout_hw]
            prof = optics.azimuthalAverage(im_cutout, returnradii=True,
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
        return im_slit

    def simulate_image(self, x, wave, blaze, matrices, im_slit, spectrum=None,
                       n_x=0, xshift=0.0, yshift=0.0, rv=0.0):
        """Simulate a spectrum on the CCD.

        Parameters
        ----------
        x,w,b,matrices: float arrays
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
        yshift: float
            NOT IMPLEMENTED
        rv: float
            Radial velocity in m/s.
        """
        # If no input spectrum, use the sun.
        if (spectrum is None) or len(spectrum) == 0:
            data = pf.getdata(os.path.join(LOCAL_DIR, 'data/ardata.fits.gz'))
            spectrum = np.array([np.append(0.35, data['WAVELENGTH']) / 1e4,
                                 np.append(0.1, data['SOLARFLUX'])])
        n_orders = x.shape[0]
        n_y = x.shape[1]
        if n_x == 0:
            n_x = n_y
        image = np.zeros((n_y, n_x))
        # Simulate the slit image within a small cutout region.
        cutout_xy = np.meshgrid(np.arange(81)-40, np.arange(7)-3)
        # Loop over orders
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
                cutout_shifted = (cutout_xy[0].copy() + int(the_x) + n_x/2,
                                  cutout_xy[1].copy() + j)
                w_ix = np.where((cutout_shifted[0] >= 0) *
                                (cutout_shifted[1] >= 0) *
                                (cutout_shifted[0] < n_x) *
                                (cutout_shifted[1] < n_y))
                cutout_shifted = (cutout_shifted[0][w_ix],
                                  cutout_shifted[1][w_ix])
                flux = np.interp(wave[i, j]*(1 + rv/299792458.0),
                                 spectrum[0], spectrum[1],
                                 left=0, right=0)
                # Rounded to the nearest microns_pix, find the co-ordinate
                # in the simulated slit image corresponding to each pixel.
                # The co-ordinate order in the matrix is (x,y).
                xy_scaled = np.dot(
                    matrices[i, j],
                    np.array([cutout_xy[0][w_ix] + int(the_x) - the_x,
                              cutout_xy[1][w_ix]]) / self.microns_pix)\
                              .astype(int)
                slit_y = xy_scaled[1] + im_slit.shape[0]/2
                slit_x = xy_scaled[0] + im_slit.shape[1]/2
                w_ix = np.where((slit_x >= 0) * (slit_y >= 0) *
                                (slit_x < im_slit.shape[1]) *
                                (slit_y < im_slit.shape[0]))
                cutout_shifted = (cutout_shifted[0][w_ix],
                                  cutout_shifted[1][w_ix])
                slit_x = slit_x[w_ix]
                slit_y = slit_y[w_ix]
                image[cutout_shifted[1], cutout_shifted[0]] += blaze[i, j] * \
                    flux * im_slit[slit_y, slit_x]
            print('Done order: {0}'.format(i + self.order_min))
        return image

    def thar_spectrum(self):
        """Calculates a ThAr spectrum.

        Returns
        -------
        wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
        """

        thar = np.loadtxt(
            os.path.join(LOCAL_DIR, 'data/mnras0378-0221-SD1.txt'),
            usecols=[0, 1, 2])
        #Create a fixed wavelength scale evenly spaced in log.
        thar_wave = 3600 * np.exp(np.arange(5e5)/5e5)
        thar_flux = np.zeros(5e5)
        #NB This is *not* perfect: we just find the nearest index of the
        #wavelength scale corresponding to each Th/Ar line.
        wave_ix = (np.log(thar[:, 1]/3600) * 5e5).astype(int)
        wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5-1).astype(int)
        thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))
        thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                                mode='same')
        # Make the peak flux equal to 10
        thar_flux /= 0.1*np.max(thar_flux)
        return np.array([thar_wave/1e4, thar_flux])

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

        # FIXME this sky spectrum isn't high enough resolution
        bgdata = np.loadtxt(
            os.path.join(LOCAL_DIR, 'data/skybg_50_10.dat'), skiprows=14)
        # Convert nm to um
        bgdata[:, 0] /= 1000

        # Calculate area per fiber
        # FIXME this is just a rough guess - should get the real numbers
        if mode == 'high':
            area = math.pi * 0.125**2
        elif mode == 'std':
            area = math.pi * 0.2**2

        # Convert phot/s/nm/arcsec^2/m^2 into phot/s/nm
        bgdata[:, 1] *= math.pi * 4**2 * area

        # Wavelength of spectrum
        bgwave = np.linspace(bgdata[:, 0].min(), bgdata[:, 0].max(), 100000)

        # Size of each step in nm
        bgstep = 1000 * (bgwave[1] - bgwave[0])

        # Convert to phot/s
        bgdata[:, 1] *= bgstep

        # Interpolate onto a regular grid
        bgflux = np.interp(bgwave, bgdata[:, 0], bgdata[:, 1])

        # Return flux and wavelength
        return np.array([bgwave, bgflux])

    def simulate_frame(self, duration=0.0, output_prefix='test_',
                       spectrum=None, xshift=0.0, yshift=0.0, rv=0.0,
                       rv_thar=0.0, flux=1e2, rnoise=3.0, gain=[1.0],
                       bias_level=10, overscan=32, namps=[1, 1],
                       use_thar=True, mode='high', add_cosmic=True,
                       add_sky=True, return_image=False, thar_flatlamp=False):
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

        spectrum: array of shape (2,n) where n is the number of
            spectral samples input spectrum.
            spectrum[0] is the wavelength array, spectrum[1] is the flux array.

        xshift: float (optional)
            x-direction (along-slit) shift.

        yshift: float (optional)
            y-direction (spectral direction) shift.

        rv: float (optional)
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
            Gain in electrons per ADU, per amplifier (or a single scalar for
            all amps).

        bias_level: float (optional)
            Bias level in electrons, per amplifier (or a single scalar for
            all amps).

        overscan: int (optional)
            number of columns of overscan

        namps: int[2] (optional)
            number of readout amps in each direction

        use_thar: bool (optional)
            Is the Thorium/Argon lamp in use?

        mode: string (optional)
            Can be 'high' or 'std' for the resolution mode.

        add_cosmic: bool (optional)
            Are cosmic rays added to the frame?

        add_sky: bool (optional)
            Is the sky background added to the frame?

        return_image: bool (optional)
            Do we return an image as an array? The fits file is always written.
        """

        x, wave, blaze, matrices = self.spectral_format_with_matrix()

        # Deal with the values that can be scalars or arrays.
        # Turn them into arrays of the right size if they're scalars.
        if np.isscalar(gain):
            gain = np.ones(namps[0]*namps[1]) * gain

        if np.isscalar(bias_level):
            bias_level = np.ones(namps[0]*namps[1]) * bias_level

        if np.isscalar(rnoise):
            rnoise = np.ones(namps[0]*namps[1]) * rnoise

        if mode == 'high':
            slit_fluxes = np.ones(19)*0.37
            slit_fluxes[6:13] = 0.78
            slit_fluxes[9] = 1.0
            slit_fluxes /= np.mean(slit_fluxes)
            im_slit = self.make_lenslets(fluxes=slit_fluxes, mode='high',
                                         llet_offset=2)
            image = self.simulate_image(x, wave, blaze, matrices, im_slit,
                                        spectrum=spectrum, n_x=self.szx,
                                        xshift=xshift, rv=rv)

            if use_thar:
                # Create an appropriately convolved Thorium-Argon spectrum
                thar_spect = self.thar_spectrum()
                '''
                plt.clf()
                plt.plot(thar_spect[0], thar_spect[1])
                plt.show()
                '''
                if thar_flatlamp:
                    thar_spect[1][:] = 10
                # Now that we have our spectrum, create the Th/Ar image.
                slit_fluxes = np.ones(1)
                im_slit2 = self.make_lenslets(fluxes=slit_fluxes, mode='high',
                                              llet_offset=0)
                image += self.simulate_image(x, wave, blaze, matrices, im_slit2,
                                             spectrum=thar_spect, n_x=self.szx,
                                             xshift=xshift, rv=rv_thar)
        else:
            print("ERROR: unknown mode.")
            raise UserWarning

        # Prevent any interpolation errors (negative flux)
        # prior to adding noise.
        image = np.maximum(image, 0)

        # Scale to photons
        image = duration * np.random.poisson(flux * image)

        if add_sky:
            # Add sky spectrum here
            sky_spect = self.sky_background(mode)
            # Put it into all the fibers equally
            if mode == 'high':
                slit_fluxes = np.ones(26)
                im_slit2 = self.make_lenslets(fluxes=slit_fluxes, mode=mode,
                                              llet_offset=0)
            elif mode == 'std':
                slit_fluxes = np.ones(10)
                # FIXME llet_offset for standard mode?
                im_slit2 = self.make_lenslets(fluxes=slit_fluxes, mode=mode,
                                              llet_offset=0)
            sky_image = self.simulate_image(x, wave, blaze, matrices, im_slit2,
                                            spectrum=sky_spect, n_x=self.szx,
                                            xshift=xshift, rv=0.0)
            sky_image = np.maximum(0, sky_image)
            image += duration * np.random.poisson(flux * sky_image)

        # Probably ignore atmospheric transmission for now

        # Scale from photons to electrons by multiplying by QE
        # FIXME do we want to simulate wavelength dependent QE?
        image *= 0.8

        if add_cosmic:
            # Add cosmic rays
            # We hard code 2.0 rays/s/cm^2
            # A shield that only allows rays that originate within 10 degrees
            # of the zenith.
            # Pixel size of 15 x 15 x 16 um
            image += cosmic.cosmic(image.shape, duration, 10, 2.0, False,
                                   [15, 15, 16])

        # Add dark current (3 e/pix/hour)
        image += np.random.poisson(np.ones_like(image) * 3.0 * duration/3600.0)

        # FIXME Multiply by flatfield map (i.e. simulate pixel-pixel
        # non-uniformity).
        # FIXME Apply non-linearity

        # For conventional axes transpose the image
        image = image.T

        # Split the image into sections for each readout amplifier
        images = np.array_split(image, namps[0])
        current_size = 0
        ccdsecx = np.zeros((namps[1], 2))
        ccdsecy = np.zeros((namps[0], 2))
        newimages = []
        for i, im_amp in enumerate(images):
            ccdsecy[i, 0] = current_size
            current_size += im_amp.shape[0]
            ccdsecy[i, 1] = current_size
            newimages.extend(np.array_split(image, namps[1], axis=1))
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

        # Add read noise for each amplifier
        images = [i+r*np.random.normal(size=i.shape)
                  for i, r in zip(images, rnoise)]

        # Divide by the gain in e/ADU for each amplifier
        images = [(a/g + b) for (a, g, b) in zip(images, gain, bias_level)]

        # Add overscan for each amplifier
        if overscan > 0:
            newimages = []
            for (i, g, b, r) in zip(images, gain, bias_level, rnoise):
                o = r * np.random.normal(size=(i.shape[0], overscan)) + b
                o /= g
                i = np.hstack((i, o))
                newimages.append(i)
            images = newimages

        # FIXME should convert to unsigned short here,
        # and deal with saturation.

        # Now create our fits image!
        # By adding the DETSIZE and DETSEC keywords we can open the
        # raw image in ds9 as an IRAF mosaic.
        hdr = pf.Header()
        hdr['DETSIZE'] = "[1:%d,1:%d]" % (image.shape[1], image.shape[0])
        hdulist = pf.HDUList(pf.PrimaryHDU(header=hdr))
        for i, im_amp in enumerate(images):
            hdr = pf.Header()
            hdr['DETSIZE'] = "[1:%d,1:%d]" % (image.shape[1], image.shape[0])
            hdr['DETSEC'] = "[%d:%d,%d:%d]" % (cxl[i]+1, cxh[i],
                                               cyl[i]+1, cyh[i])
            hdr['DATASEC'] = "[%d:%d,%d:%d]" % (1, im_amp.shape[1]-overscan,
                                                1, im_amp.shape[0])
            hdr['TRIMSEC'] = "[%d:%d,%d:%d]" % (1, im_amp.shape[1]-overscan,
                                                1, im_amp.shape[0])
            if overscan > 0:
                hdr['BIASSEC'] = "[%d:%d,%d:%d]" % (im_amp.shape[1]-overscan+1,
                                                    im_amp.shape[1], 1, im_amp.shape[0])
            hdulist.append(pf.ImageHDU(data=im_amp, header=hdr))
        hdulist.writeto(output_prefix + self.arm + '.fits', clobber=True)

        if return_image:
            return images