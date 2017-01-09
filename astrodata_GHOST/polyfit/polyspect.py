"""This is a simple simulation and extraction definition code for RHEA.
The key is to use Tobias's spectral fitting parameters and define
the same kinds of extraction arrays as pyghost.
"""
from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider, Button, RadioButtons
import os
import scipy.optimize as op
import pdb



class Polyspect(object):
    """A class containing tools common for any spectrograph.
    This class should be inhereted by the specific spectrograph module.
    Contains functions related to polynomial modelling of spectrograph orders.
    """

    def __init__(self,mode,m_ref,szx,szy,m_min,m_max,
                 transpose,extra_rot,nl,fiber_separation,profile_sigma):
        """ Initialisation function for this class"""

        # All necessary parameters are listed here and initialised by the
        # parent class
        self.mode=mode
        self.m_ref=m_ref
        self.szx=szx
        self.szy=szy
        self.m_min=m_min
        self.m_max=m_max
        self.transpose=transpose
        self.extra_rot=extra_rot
        self.nl=nl
        self.fiber_separation=fiber_separation
        self.profile_sigma=profile_sigma
               
    
    def wave_fit_resid(self, params, ms, waves, ys, ydeg=3, xdeg=3):
        """A fit function for read_lines_and_fit (see that function for details)
        to be used in scipy.optimize.leastsq as the minimisation function. 
        The same function is used in fit_to_x, but in that case 
        "waves" is replaced by "xs".

        Parameters
        ----------
        
        params: float array
            2D array containing the polynomial coefficients that will form the 
            model to be compared with the real data. 
        ms: int array
            The order numbers for the residual fit repeated ys times.
        waves: float array
            This parameter is the argument for the second set of polynomials. 
            Either the wavelengths or the x positions on the chip
        ys: float array
            This is an ms x y sized array with pixel indeces on y direction
            for each order.
        ydeg: int
            Polynomial degree as a function of order
        xdeg: int
            Polynomial degree as a function of y

        Returns
        -------

        The residual between the model and data supplied.
        """
        if np.prod(params.shape) != (xdeg + 1) * (ydeg + 1):
            print("Parameters are flattened - xdeg and ydeg must be correct!")
            raise UserWarning
        params = params.reshape((ydeg + 1, xdeg + 1))
        if (len(ms) != len(ys)):
            print("ms and ys must all be the same length!")
            raise UserWarning
        mp = self.m_ref / ms - 1

        ps = np.empty((len(ms), ydeg + 1))
        # Find the polynomial coefficients for each order.
        for i in range(ydeg + 1):
            polyq = np.poly1d(params[i, :])
            ps[:, i] = polyq(mp)
        wave_mod = np.empty(len(ms))
        for i in range(len(ms)):
            polyp = np.poly1d(ps[i, :])
            wave_mod[i] = polyp(ys[i] - self.szy // 2)
        return wave_mod - waves

    def read_lines_and_fit(self, init_mod, pixdir='',
                           outdir='./', ydeg=3, xdeg=3, residfile='resid.txt'):
        """Read in a series of text files that have a (Wavelength, pixel)
        format in file names like order99.txt and order100.txt.
        Fit an nth order polynomial to the wavelength as a function
        of pixel value. THIS IS A PLACE HOLDER FUNCTION NOT PART OF THE
        LATEST RELEASE. THIS WILL BE USED FOR WAVELENGTH FITTING.

        The functional form is:

        ..math::
        
            wave = p_0(m) + p_1(m)*y' + p_2(m)*y'^2 + ...)

        with :math: `y^{prime} = y - y_{\rm middle}`, and:

        .. math::

            p_0(m) = q_{00} + q_{01} * m' + q_{02} * m'^2 + ...

        with :math: `m^{prime} = m_{\rm ref}/m - 1`

        This means that the simplest spectrograph model should have:
        :math: `q_{00}` : central wavelength or order m_ref
        :math: `q_{01}` : central wavelength or order m_ref
        :math: `q_{10}` : central_wavelength/R_pix, with R_pix the resolving power / pixel.
        :math: `q_{11}` : central_wavelength/R_pix, with R_pix the resolving power / pixel.
        ... with everything else approximately zero.

        Parameters
        ----------
        
        init_mod: 2D array
            initial model parameters.
        pixdir: string
            Alternative location for the directory containing the arclines info
            THIS IS CURRENTLY NOT COMPLIENT WITH THE PROCEDURE WE WANT. ARC 
            LINES MUST BE FED WITHOUT A KNOWLEDGE OF THE FILE LOCATION.
        xdeg, ydeg: int
            Order of polynomial
        residfile: string
            Residual file output name

        Returns
        -------
        params: float array
            Fitted parameters
        wave_and_resid: float array
            Wavelength and fit residuals.
        """
        if (len(pixdir) == 0):
            pixdir = self.model_location
        else:
            if type(pixdir) is not str:
                return 'Location for arcline files invalid'
        # The next loop reads in wavelengths from a file.
        # To make this neater, it could be a function that overrides this
        # base class.
        # This code needs more work since we will only have one format
        # for a GHOST arc line list
        if os.path.exists(pixdir + "arclines.txt"):
            lines = np.loadtxt(pixdir + "arclines.txt")
            ms = lines[:, 3]
            waves = lines[:, 0]
            ys = lines[:, 1]
        else:
            ms = np.array([])
            waves = np.array([])
            ys = np.array([])
            for m in range(self.m_min, self.m_max + 1):
                fname = pixdir + "order{0:d}.txt".format(m)
                try:
                    pix = np.loadtxt(fname)
                except:
                    print("Error: arc line files don't exist!")
                    raise UserWarning
                ms = np.append(ms, m * np.ones(pix.shape[0]))
                waves = np.append(waves, pix[:, 0])
                # Tobias's definition of the y-axis is backwards compared to
                # python.
                ys = np.append(ys, self.szy - pix[:, 1])

        init_resid = self.wave_fit_resid(init_mod, ms, waves, ys, ydeg=ydeg,
                                         xdeg=xdeg)
        bestp = op.leastsq(self.wave_fit_resid, init_mod, args=(ms, waves, ys,
                                                               ydeg, xdeg))
        final_resid = self.wave_fit_resid(bestp[0], ms, waves, ys, ydeg=ydeg,
                                          xdeg=xdeg)
        # Output the fit residuals.
        wave_and_resid = np.array([waves, ms, final_resid]).T
        print("Fit residual RMS (Angstroms): {0:6.3f}".format(
            np.std(final_resid)))
        params = bestp[0].reshape((ydeg + 1, xdeg + 1))
        return params, wave_and_resid

    
    def spectral_format(self, xoff=0.0, yoff=0.0, ccd_centre={}, wparams=None,
                        xparams=None, img=None):
        """Create a spectrum, with wavelengths sampled in 2 orders based on
           a pre-existing wavelength and x position polynomial model.
           This code takes the polynomial model and calculates the result as 
           a function of order number scaled to the reference order and then
           as a function of y position. 
           Optionally a file can be supplied for the model to be overlayed 
           on top of.

        Parameters
        ----------
        xoff: float
            An input offset from the field center in the slit plane in
            mm in the x (spatial) direction.
        yoff: float
            An input offset from the field center in the slit plane in
            mm in the y (spectral) direction.
        ccd_centre: dict
            An input describing internal parameters for the angle of
            the center of the CCD. To run this program multiple times
            with the same co-ordinate system, take the returned
            ccd_centre and use it as an input.
        wparams: float array (optional)
            2D array with polynomial parameters for wavelength scale
        xparams: float array (optional)
            2D array with polynomial parameters for x scale
        img: 2D array (optional)
            2D array containing an image. This function
            uses this image and over plots the created position model.

        Returns
        -------
        x:  (nm, ny) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel and
            each order (m).
        wave: (nm, ny) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (nm, ny) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        ccd_centre: dict
            NOT YET IMPLEMENTED
            Parameters of the internal co-ordinate system describing the
            center of the CCD.
        """

        # Now lets interpolate onto a pixel grid rather than the arbitrary
        # wavelength grid we began with.
        nm = self.m_max - self.m_min + 1
        x_int = np.zeros((nm, self.szy))
        wave_int = np.zeros((nm, self.szy))
        blaze_int = np.zeros((nm, self.szy))

        if (xparams is None) and (wparams is None):
            return 'Must provide at least one of xparams or wparams'
        ys = np.arange(self.szy)
        # Loop through m
        for m in np.arange(self.m_min, self.m_max + 1):
            # First, sort out the wavelengths
            mp = self.m_ref / m - 1

            if wparams is not None:
                # Find the polynomial coefficients for each order.
                ydeg = wparams.shape[0] - 1
                ps = np.empty((ydeg + 1))
                for i in range(ydeg + 1):
                    polyq = np.poly1d(wparams[i, :])
                    ps[i] = polyq(mp)
                polyp = np.poly1d(ps)
                wave_int[m - self.m_min, :] = polyp(ys - self.szy // 2)

            if xparams is not None:
                # Find the polynomial coefficients for each order.
                ydeg = xparams.shape[0] - 1
                ps = np.empty((ydeg + 1))
                for i in range(ydeg + 1):
                    polyq = np.poly1d(xparams[i, :])
                    ps[i] = polyq(mp)
                polyp = np.poly1d(ps)
                x_int[m - self.m_min, :] = polyp(ys - self.szy // 2)

            # Finally, the blaze
            wcen = wave_int[m - self.m_min, self.szy / 2]
            disp = wave_int[m - self.m_min, self.szy / 2 + 1] - wcen
            order_width = (wcen / m) / disp
            blaze_int[m - self.m_min, :] = np.sinc((ys - self.szy / 2) /
                                                   order_width)**2

        # Plot this if we have an image file 
        if (img is not None) and (xparams is not None):
            if not self.transpose:
                img = img.T
            plt.clf()
            plt.imshow(np.arcsinh((img - np.median(img)) / 100), aspect='auto',
                       interpolation='nearest', cmap=cm.gray)
            plt.axis([0, img.shape[1], img.shape[0], 0])
            plt.plot(x_int.T + + self.szx // 2)

        return x_int, wave_int, blaze_int

    def adjust_x(self, old_x, image, num_xcorr=21):
        """Adjust the x pixel value based on an image and an initial
            array from spectral_format().
            This function performs a cross correlation between a pixel map
            and a given 2D array image and calculates a shift along the spatial
            direction. The result is used to inform an initial global shift for
            the fitting procedure.
            Only really designed for a single fiber flat, single fiber science
            image or a convolution map.
            This is a helper routine for fit_x_to_image.

        Parameters
        ----------
        old_x: numpy array
            An old x pixel array
        image: numpy array
            A 2D image array to be used as the basis for the adjustment.

        Returns
        -------
        A new value of the x array.
        """
        # Create an array with a single pixel with the value 1.0 at the
        # expected peak of each order.
        single_pix_orders = np.zeros(image.shape)
        xy = np.meshgrid(np.arange(old_x.shape[0]), np.arange(old_x.shape[1]))
        single_pix_orders[np.round(xy[1]).astype(int),
                          np.round(old_x.T + self.szx // 2).astype(int)] = 1.0

        # Make an array of cross-correlation values.
        xcorr = np.zeros(num_xcorr)
        for i in range(num_xcorr):
            xcorr[i] = np.sum(np.roll(single_pix_orders, i - num_xcorr // 2,
                                      axis=1) * image)

        # Based on the maximum cross-correlation, adjust the model x values.
        the_shift = np.argmax(xcorr) - num_xcorr // 2

        return old_x + the_shift

    def fit_x_to_image(self, data, xparams, decrease_dim=10, search_pix=20,
                       xdeg=4,inspect=False,clobber=False):
        """Fit a "tramline" map. Note that an initial map has to be pretty
        close, i.e. within "search_pix" everywhere. To get within search_pix
        everywhere, a simple model with a few paramers is fitted manually.
        This can be done with the GUI using the adjust_model function and 
        finally adjusted with the adjust_x function

        Parameters
        ----------
        data: numpy array
            The image of a single reference fiber to fit to. Typically 
            the result of the convolution.
        xparams: float array
            The polynomial parameters to be fitted. 
        decrease_dim: int
            Median filter by this amount in the dispersion direction and
            decrease the dimensionality of the problem accordingly.
            This helps with both speed and robustness.
        search_pix: int
            Search within this many pixels of the initial model.
        xdeg: int
            Polynomial degree. This parameter is probably not needed.
        inspect: bool
            If true, once fit is done the adjust_model function
            is called so that the user can inspect the result of the
            fit and decide if it is good enough. 
        clobber: bool
            If true this will automatically replace the existing model
            with the fitted one.
        """
        xx, wave, blaze = self.spectral_format(xparams=xparams)
        if self.transpose:
            image=data.T
        else:
            image=data
        # Now use the adjust function to figure out a global shift in
        # the spatial direction
        xs = self.adjust_x(xx, image)
        
        if image.shape[0]%decrease_dim!=0:
            return "Can not decrease image dimention by this amount. "+\
            "Please check if the image size in the spectral dimention is "+\
            "exactly divisible by this amount."
        # Median-filter in the dispersion direction.
        # This process will 'collapse' the image in the spectral direction
        # and make sure the fit is faster. 
        image_med = image.reshape((image.shape[0] // decrease_dim,
                                   decrease_dim, image.shape[1]))
        image_med = np.median(image_med, axis=1)
        my = np.meshgrid(np.arange(xx.shape[1]), np.arange(xx.shape[0]) +
                         self.m_min)
        ys = my[0]
        ys = np.average(ys.reshape(xs.shape[0], xs.shape[1] // decrease_dim,
                                   decrease_dim), axis=2)
        xs = np.average(xs.reshape(xs.shape[0], xs.shape[1] // decrease_dim,
                                   decrease_dim), axis=2)

        # Now go through and find the peak pixel values.
        # Do this by searching for the maximum value along the
        # order for search_pix on either side of the initial
        # model pixels in the spatial direction. 
        for i in range(xs.shape[0]):  # Go through each order...
            for j in range(xs.shape[1]):
                xi = int(np.round(xs[i, j]))
                peakpix = image_med[j, self.szx // 2 + xi -
                                    search_pix:self.szx // 2 +
                                    xi + search_pix + 1]
                xs[i, j] += np.argmax(peakpix) - search_pix

        fitted_params=self.fit_to_x(xs, xparams, ys=ys, xdeg=xdeg)
        if inspect:
            #This will plot the result of the fit once successful so
            #the user can inspect the result.
            plt.imshow((data - np.median(data)) / 1e2)
            x_int,wave_int,blaze_int = self.spectral_format(wparams=None,
                                                            xparams=fitted_params)
            y = np.meshgrid(np.arange(data.shape[1]),
                            np.arange(x_int.shape[0]))[0]
            plt.plot(y,x_int + data.shape[0]//2,
                     color='green', linestyle='None', marker='.')
            plt.show()

        return fitted_params

    def fit_to_x(self, x_to_fit, init_mod, ydeg=2,
                 xdeg=4, ys=[], decrease_dim=1):
        """Fit to an (nm,ny) array of x-values.

        The functional form is:
        
        .. math::
        
            x = p_0(m) + p_1(m)*y' + p_2(m)*y'^2 + ...)

        with :math:`y^{prime} = y - y_{\rm middle}`, and:
        
        .. math::
        
            p_0(m) = q_{00} + q_{01} * m' + q_{02} * m'^2 + ...

        with :math:`mp = m_{\rm ref}/m - 1

        This means that the simplest spectrograph model should have:
        :math:`q_{00}` : central order y pixel
        :math:`q_{01}`:  spacing between orders divided by the number of orders
        ... with everything else approximately zero.

        Parameters
        ----------
        x_to_fit: float array
            x values to fit. This should be an (orders,y) shape array.
        init_mod_file: float array
            Initial model parameters
        xdeg, ydeg: int
            Order of polynomial
        y: float array
            Y positions on the CCD
        decrease_dim: int
            The factor of decreased dimentionality for the fit. 
            This needs to be an exact factor of the y size.        
        """

        # Create an array of y and m values.
        xs = x_to_fit.copy()
        my = np.meshgrid(np.arange(xs.shape[1]), np.arange(xs.shape[0]) +
                         self.m_min)
        if (len(ys) == 0):
            ys = my[0]
        ms = my[1]

        # Allow a dimensional decrease, for speed
        if (decrease_dim > 1):
            ms = np.average(ms.reshape(xs.shape[0], xs.shape[
                            1] // decrease_dim, decrease_dim), axis=2)
            ys = np.average(ys.reshape(xs.shape[0], xs.shape[1] //
                                       decrease_dim, decrease_dim), axis=2)
            xs = np.average(xs.reshape(xs.shape[0], xs.shape[1] //
                                       decrease_dim, decrease_dim), axis=2)

        # Flatten arrays
        ms = ms.flatten()
        ys = ys.flatten()
        xs = xs.flatten()

        # Do the fit!
        print("Fitting (this can sometimes take a while...)")
        init_resid = self.wave_fit_resid(
            init_mod, ms, xs, ys, ydeg=ydeg, xdeg=xdeg)
        bestp = op.leastsq(self.wave_fit_resid, init_mod,
                           args=(ms, xs, ys, ydeg, xdeg))
        final_resid = self.wave_fit_resid(
            bestp[0], ms, xs, ys, ydeg=ydeg, xdeg=xdeg)
        params = bestp[0].reshape((ydeg + 1, xdeg + 1))
        print(init_resid, final_resid)

        return params

    def spectral_format_with_matrix(self):
        """Create a spectral format, including a detector to slit matrix at
           every point.

        Returns
        -------
        x: (nm, ny) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel
            and each order (m).
        w: (nm, ny) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (nm, ny) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        matrices: (nm, ny, 2, 2) float array
            2x2 slit rotation matrices, mapping output co-ordinates back
            to the slit.
        """
        x, w, b = self.spectral_format()
        matrices = np.zeros((x.shape[0], x.shape[1], 2, 2))
        amat = np.zeros((2, 2))

        for i in range(x.shape[0]):  # i is the order
            for j in range(x.shape[1]):
                # Create a matrix where we map input angles to output
                # coordinates.
                slit_microns_per_det_pix = self.slit_microns_per_det_pix_first + \
                    float(i) / x.shape[0] * (self.slit_microns_per_det_pix_last - \
                                             self.slit_microns_per_det_pix_first)
                amat[0, 0] = 1.0 / slit_microns_per_det_pix
                amat[0, 1] = 0
                amat[1, 0] = 0
                amat[1, 1] = 1.0 / slit_microns_per_det_pix
                # Apply an additional rotation matrix. If the simulation was
                # complete, this wouldn't be required.
                r_rad = np.radians(self.extra_rot)
                dy_frac = (j - x.shape[1] / 2.0) / (x.shape[1] / 2.0)
                extra_rot_mat = np.array([[np.cos(r_rad * dy_frac),
                                           np.sin(r_rad * dy_frac)],
                                          [-np.sin(r_rad * dy_frac),
                                           np.cos(r_rad * dy_frac)]])
                amat = np.dot(extra_rot_mat, amat)
                # We actually want the inverse of this (mapping output
                # coordinates back onto the slit.
                matrices[i, j, :, :] = np.linalg.inv(amat)
        return x, w, b, matrices

    def slit_flat_convolve(self, flat, slit_profile=None):
        """Function that takes a flat field image and a slit profile and
           convolves the two in 2D. Returns result of convolution, which
           should be used for tramline fitting. 

        Parameters
        ----------
        flat: float array
            A flat field image from the spectrograph
        slit_profile: float array
            A 1D array with the slit profile fiber amplitudes. If none is
            supplied this function will assume identical fibers and create
            one to be used in the convolution based on default parameters 
            specified in the ghost class.

        Returns
        -------
        flat_conv: float array
            Returns the convolved 2D array.

        """
        if not slit_profile:
            # At this point create a slit profile
            # Create a x baseline for convolution and assume a FWHM for profile
            x = flat.shape[0]
            profilex = np.arange(x) - x // 2
            # Now create a model of the slit profile
            mod_slit = np.zeros(x)
            if self.mode == 'high':
                nfibers = 26
            else:
                nfibers = self.nl

            for i in range(-nfibers // 2, nfibers // 2):
                mod_slit += np.exp(-(profilex - i * self.fiber_separation)**2 /
                                   2.0 / self.profile_sigma**2)
        else:
            mod_slit=slit_profile

        # Normalise the slit model and fourier transform for convolution
        mod_slit /= np.sum(mod_slit)
        mod_slit_ft = np.fft.rfft(np.fft.fftshift(mod_slit))
        # Fourier the flat for convolution
        im_fft = np.fft.rfft(flat, axis=0)
        flat_conv = np.zeros_like(im_fft)
        # Now convolved in 2D
        for i in range(im_fft.shape[1]):
            flat_conv[:, i] = im_fft[:, i] * mod_slit_ft
        flat_conv = np.fft.irfft(flat_conv, axis=0)
        return flat_conv

    def adjust_model(self, data, wparams=None,
                      xparams=None,convolve=True,
                      percentage_variation=10,vary_wrt_max=True):
        """Function that uses matplotlib slider widgets to adjust a polynomial
        model overlaid on top of a flat field image. In practice this will be
        overlaid on top of the result of convolving the flat with a slit
        profile in 2D which just reveals the location of the middle of the
        orders.

        Parameters
        ----------
        data: float array
            an array containing data to be used as a visual comparison of the
            model
        wparams: float array (optional)
            2D array containing the initial wavelength
            model parameters.
        xparams: float array (optional)
            2D array containing the initial order location model parameters.
        convolve: boolean
            A boolean indicating whether the data provided should be convolved 
            with a slit profile model. True if data provided is a multi fiber 
            flatfield.
        percentage_variation: int
            How much should the percentage adjustment in each bin as a function
            of the parameter. Default is 10%
        vary_wrt_max: bool
            Vary all parameters intelligently with a scaling of the maximum variation, 
            rather than just a percentage of each.

        Returns
        -------
        Either a success or an error message.

        """
        #Must provide at least one of the model parameters
        if (xparams is None) and (wparams is None):
            return 'Must provide at least one of xparams or wparams'
        
        nx = data.shape[0]
        # Start by setting up the graphical part
        fig, ax = plt.subplots()
        axcolor = 'lightgoldenrodyellow'
        # Grab the model to be plotted
        x_int, wave_int, blaze_int = self.spectral_format(wparams=wparams,
                                                           xparams=xparams)
        y = np.meshgrid(np.arange(data.shape[1]), np.arange(x_int.shape[0]))[0]
        # The data must be flattened for the sliders to work.
        # Then plot it! 
        l, = plt.plot(y.flatten()[::10], x_int.flatten()[::10] + nx // 2,
                      color='green', linestyle='None', marker='.')

        if convolve:
            data = self.slit_flat_convolve(flat=data)
        # Now over plot the image.
        ax.imshow((data - np.median(data)) / 1e2)

        # Create a second window for cliders.
        slide_fig = plt.figure()

        # This function is executed on each slider change.
        # spectral_format is updated.
        def update(val):
            for i in range(npolys):
                for j in range(polyorder):
                    xparams[i, j] = sliders[i][j].val
            x, wave, blaze = self.spectral_format(wparams=wparams,
                                                   xparams=xparams)
            l.set_ydata(x.flatten()[::10] + nx // 2)
            fig.canvas.draw_idle()

        polyorder = xparams.shape[1]
        npolys = xparams.shape[0]
        # Now we start putting sliders in depending on number of parameters
        height = 1. / (npolys * 2)
        width = 1. / (polyorder * 2)
        # Use this to adjust in a percentage how much to let each parameter
        # vary
        frac_xparams = np.absolute(xparams * (percentage_variation / 100))
        if vary_wrt_max:
            for i in range(npolys):
                frac_xparams[i] = np.max(frac_xparams[-1])/(nx/2.0)**(npolys-1-i)
        axq = [[0 for x in range(polyorder)] for y in range(npolys)]
        sliders = [[0 for x in range(polyorder)] for y in range(npolys)]
        # Now put all the sliders in the new figure based on position in the
        # array
        for i in range(npolys):
            for j in range(polyorder):
                left = j * width * 2
                bottom = 1 - (i + 1) * height * 2 + height
                axq[i][j] = plt.axes([left, bottom, width, height],
                                     axisbg=axcolor)
                if xparams[i, j] == 0:
                    sliders[i][j] = Slider(axq[i][j],
                                           'test' + str(i) + str(j), 0, 0.1,
                                           valinit=xparams[i, j])
                else:
                    sliders[i][j] = Slider(axq[i][j], 'test' + str(i) + str(j),
                                           xparams[i, j] - frac_xparams[i, j],
                                           xparams[i, j] + frac_xparams[i, j],
                                           valinit=xparams[i, j])
                plt.legend(loc=3)
                sliders[i][j].on_changed(update)

        # Have a button to output the current model to the correct file
        submitax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(submitax, 'Submit', color=axcolor, hovercolor='0.975')

        # This is the triggered function on submit.
        # Currently only works for the xmod but should be generalised
        def submit(event):
            plt.close('all')
            return xparams
            
        button.on_clicked(submit)

        plt.show()
        return xparams
