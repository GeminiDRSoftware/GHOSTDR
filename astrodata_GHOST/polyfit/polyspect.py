"""This is a simple simulation and extraction definition code for RHEA.
The key is to use Tobias's spectral fitting parameters and define
the same kinds of extraction arrays as pyghost.
"""
from __future__ import division, print_function

import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider, Button
import scipy.optimize as op

# pylint: disable=maybe-no-member, too-many-instance-attributes

class Polyspect(object):
    """A class containing tools common for any spectrograph.
    This class should be inhereted by the specific spectrograph module.
    Contains functions related to polynomial modelling of spectrograph orders.
    """

    def __init__(self, m_ref, szx, szy, m_min, m_max, transpose):
        """ Initialisation function for this class"""

        # All necessary parameters are listed here and initialised by the
        # parent class
        self.m_ref = m_ref
        self.szx = szx
        self.szy = szy
        self.m_min = m_min
        self.m_max = m_max# 
        #True if the spectral dispersion dimention is over the x (column) axis
        self.transpose = transpose

    def evaluate_poly(self, params, mprime, y_values):
        """ Function used to evaluate a polynomial of polynomials at specific
        points. This is a function designed to avoid code repetition.

        Parameters
        ----------

        params: float array
            Model parameters with the coefficients to evaluate.
        mprime: float
            This is the mprime as defined in the spectral format function. This
            parameter is what the first order of polynomials are evaluated
            against
        y_values: float array
            This is the value or array of values for the second order polynomial
            functions to be aveluated against."""
        #params needs to be a np.array
        if not isinstance(params,np.ndarray):
            raise TypeError('Please provide params as a numpy float array')
        if not isinstance(mprime,float):
            raise TypeError('Please ensure mprime is a float number')
        if not isinstance(y_values,np.ndarray):
            raise TypeError('Please provide y_values as a numpy float array')
        ydeg = params.shape[0] - 1
        if params.ndim==1:
            polyp = np.poly1d(params)
            evaluation = polyp(mprime)
        else:            
            polynomials = np.empty((ydeg + 1))
            for i in range(ydeg + 1):
                polyq = np.poly1d(params[i, :])
                polynomials[i] = polyq(mprime)
            polyp = np.poly1d(polynomials)
            evaluation = polyp(y_values - self.szy // 2)
        return evaluation
    

    def wave_fit_resid(self, params, orders, waves, y_values, ydeg=3, xdeg=3):
        """A fit function for read_lines_and_fit (see that function for details)
        to be used in scipy.optimize.leastsq as the minimisation function.
        The same function is used in fit_to_x, but in that case
        "waves" is replaced by "xs".

        Parameters
        ----------

        params: float array
            2D array containing the polynomial coefficients that will form the
            model to be compared with the real data.
        orders: int array
            The order numbers for the residual fit repeated ys times.
        waves: float array
            This parameter is the argument for the second set of polynomials.
            Either the wavelengths or the x positions on the chip
        y_values: float array
            This is an orders x y sized array with pixel indeces on y direction
            for each order.
        ydeg: int
            Polynomial degree as a function of order
        xdeg: int
            Polynomial degree as a function of y

        Returns
        -------

        The residual between the model and data supplied.
        """
        #params needs to be a np.array
        if not isinstance(params,np.ndarray):
            raise TypeError('Please provide params as a numpy float array')
        if not isinstance(orders,np.ndarray):
            raise TypeError('Please ensure orders is a numpy float array')
        if not isinstance(y_values,np.ndarray):
            raise TypeError('Please provide y_values as a numpy float array')
        if not isinstance(waves,np.ndarray):
            raise TypeError('Please provide waves as a numpy float array')

        if np.prod(params.shape) != (xdeg + 1) * (ydeg + 1):
            print("Parameters are flattened - xdeg and ydeg must be correct!")
            raise UserWarning
        params = params.reshape((ydeg + 1, xdeg + 1))
        if len(orders) != len(y_values):
            print("orders and y_values must all be the same length!")
            raise UserWarning
        mprime = self.m_ref / orders - 1

        polynomials = np.empty((len(orders), ydeg + 1))
        # Find the polynomial coefficients for each order.
        for i in range(ydeg + 1):
            polyq = np.poly1d(params[i, :])
            polynomials[:, i] = polyq(mprime)
        wave_mod = np.empty(len(orders))
        for i in range(len(orders)):
            polyp = np.poly1d(polynomials[i, :])
            wave_mod[i] = polyp(y_values[i] - self.szy // 2)
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
        :math: `q_{00}` : central wavelength of order m_ref
        :math: `q_{01}` : central wavelength of order m_ref
        :math: `q_{10}` : central_wavelength/R_pix,
        with R_pix the resolving power / pixel.
        :math: `q_{11}` : central_wavelength/R_pix,
        with R_pix the resolving power / pixel.
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
        if len(pixdir) == 0:
            print('THIS REQUIRES FURTHER WORK.')
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
            orders = lines[:, 3]
            waves = lines[:, 0]
            y_values = lines[:, 1]
        else:
            orders = np.array([])
            waves = np.array([])
            y_values = np.array([])
            for order in range(self.m_min, self.m_max + 1):
                fname = pixdir + "order{0:d}.txt".format(order)
                try:
                    pix = np.loadtxt(fname)
                except:
                    print("Error: arc line files don't exist!")
                    raise UserWarning
                orders = np.append(orders, order * np.ones(pix.shape[0]))
                waves = np.append(waves, pix[:, 0])
                # Tobias's definition of the y-axis is backwards compared to
                # python.
                y_values = np.append(y_values, self.szy - pix[:, 1])

        init_resid = self.wave_fit_resid(init_mod, orders, waves,
                                         y_values, ydeg=ydeg, xdeg=xdeg)
        bestp = op.leastsq(self.wave_fit_resid, init_mod, args=(orders, waves,
                                                                y_values,
                                                                ydeg, xdeg))
        final_resid = self.wave_fit_resid(bestp[0], orders, waves, y_values,
                                          ydeg=ydeg, xdeg=xdeg)
        # Output the fit residuals.
        wave_and_resid = np.array([waves, orders, final_resid]).T
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
        norders = self.m_max - self.m_min + 1
        x_int = np.zeros((norders, self.szy))
        wave_int = np.zeros((norders, self.szy))
        blaze_int = np.zeros((norders, self.szy))

        if (xparams is None) and (wparams is None):
            raise UserWarning('Must provide at least one of xparams or wparams')
        if (xparams is not None) and (not isinstance(xparams,np.ndarray) ):
            raise UserWarning('xparams provided with invalid format')
        if (wparams is not None) and (not isinstance(wparams,np.ndarray) ):
            raise UserWarning('xparams provided with invalid format')
        y_values = np.arange(self.szy)
        # Loop through order
        for order in np.arange(self.m_min, self.m_max + 1):
            # First, sort out the wavelengths
            mprime = self.m_ref / order - 1

            if wparams is not None:
                # Find the polynomial coefficients for each order.
                wave_int[order - self.m_min, :] = self.evaluate_poly(wparams,
                                                                     mprime,
                                                                     y_values)

            if xparams is not None:
                # Find the polynomial coefficients for each order.
                x_int[order - self.m_min, :] = self.evaluate_poly(xparams,
                                                                  mprime,
                                                                  y_values)

            # Finally, the blaze
            wcen = wave_int[int(order - self.m_min), int(self.szy / 2)]
            disp = wave_int[int(order - self.m_min), int(self.szy / 2 + 1)] - wcen
            order_width = (wcen / order) / disp
            blaze_int[order - self.m_min, :] = np.sinc((y_values - self.szy / 2)
                                                       / order_width)**2


        # Plot this if we have an image file
        if (img is not None) and (xparams is not None):
            if not isinstance(img,ndarray):
                raise UserWarning('img must be numpy array')
            if img.ndim != 2:
                raise UserWarning('Image array provided is not a 2 dimensional\
                array')
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
        if not isinstance(old_x, np.ndarray):
            raise TypeError('old_x must be a numpy array')
        if not isinstance(image, np.ndarray):
            raise TypeError('image must be a numpy array')
        if image.ndim != 2:
            raise UserWarning('image array must be 2 dimentional')
        
        # Create an array with a single pixel with the value 1.0 at the
        # expected peak of each order.
        single_pix_orders = np.zeros(image.shape)
        xygrid = np.meshgrid(np.arange(old_x.shape[0]),
                             np.arange(old_x.shape[1]))
        single_pix_orders[np.round(xygrid[1]).astype(int),
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
                       xdeg=4, inspect=False):
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
        xbase, wave, blaze = self.spectral_format(xparams=xparams)
        if self.transpose:
            image = data.T
        else:
            image = data
        # Now use the adjust function to figure out a global shift in
        # the spatial direction
        x_values = self.adjust_x(xbase, image)

        if image.shape[0] % decrease_dim != 0:
            return "Can not decrease image dimention by this amount. " +\
                "Please check if the image size in the spectral dimention is " +\
                "exactly divisible by this amount."
        # Median-filter in the dispersion direction.
        # This process will 'collapse' the image in the spectral direction
        # and make sure the fit is faster.
        image_med = image.reshape((image.shape[0] // decrease_dim,
                                   decrease_dim, image.shape[1]))
        image_med = np.median(image_med, axis=1)
        order_y = np.meshgrid(np.arange(xbase.shape[1]),                        #pylint: disable=maybe-no-member
                              np.arange(xbase.shape[0]) + self.m_min)           #pylint: disable=maybe-no-member
        y_values = order_y[0]
        y_values = np.average(y_values.reshape(x_values.shape[0],
                                               x_values.shape[1] //
                                               decrease_dim,
                                               decrease_dim), axis=2)
        x_values = np.average(x_values.reshape(x_values.shape[0],
                                               x_values.shape[1] //
                                               decrease_dim,
                                               decrease_dim), axis=2)

        # Now go through and find the peak pixel values.
        # Do this by searching for the maximum value along the
        # order for search_pix on either side of the initial
        # model pixels in the spatial direction.
        for i in range(x_values.shape[0]):  # Go through each order...
            for j in range(x_values.shape[1]):                                  #pylint: disable=maybe-no-member
                xind = int(np.round(x_values[i, j]))
                peakpix = image_med[j, self.szx // 2 + xind -
                                    search_pix:self.szx // 2 +
                                    xind + search_pix + 1]
                x_values[i, j] += np.argmax(peakpix) - search_pix

        fitted_params = self.fit_to_x(x_values, xparams, y_values=y_values,
                                      xdeg=xdeg)
        if inspect:
            # This will plot the result of the fit once successful so
            # the user can inspect the result.
            plt.imshow((data - np.median(data)) / 1e2)
            x_int, wave_int, blaze_int = self.spectral_format(wparams=None,
                                                              xparams=
                                                              fitted_params)
            ygrid = np.meshgrid(np.arange(data.shape[1]),
                                np.arange(x_int.shape[0]))[0]                   #pylint: disable=maybe-no-member
            plt.plot(ygrid, x_int + data.shape[0] // 2,
                     color='green', linestyle='None', marker='.')
            plt.show()

        return fitted_params

    def fit_to_x(self, x_to_fit, init_mod, ydeg=2,
                 xdeg=4, y_values=[], decrease_dim=1):
        """Fit to an (norders,ny) array of x-values.

        The functional form is:

        .. math::

            x = p_0(m) + p_1(m)*y' + p_2(m)*y'^2 + ...)

        with :math:`y^{prime} = y - y_{\rm middle}`, and:

        .. math::

            p_0(m) = q_{00} + q_{01} * m' + q_{02} * m'^2 + ...

        with :math:`mprime = m_{\rm ref}/m - 1

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
        y_values: float array
            Y positions on the CCD
        decrease_dim: int
            The factor of decreased dimentionality for the fit.
            This needs to be an exact factor of the y size.
        """

        # Create an array of y and m values.
        x_values = x_to_fit.copy()
        order_y = np.meshgrid(np.arange(x_values.shape[1]),
                              np.arange(x_values.shape[0]) + self.m_min)
        if len(y_values) == 0:
            y_values = order_y[0]
        orders = order_y[1]

        # Allow a dimensional decrease, for speed
        if decrease_dim > 1:
            orders = np.average(orders.reshape(x_values.shape[0],
                                               x_values.shape[1] //
                                               decrease_dim, decrease_dim),
                                axis=2)
            y_values = np.average(y_values.reshape(x_values.shape[0],           #pylint: disable=maybe-no-member
                                                   x_values.shape[1] //
                                                   decrease_dim, decrease_dim),
                                  axis=2)
            x_values = np.average(x_values.reshape(x_values.shape[0],
                                                   x_values.shape[1] //
                                                   decrease_dim, decrease_dim),
                                  axis=2)

        # Flatten arrays
        orders = orders.flatten()
        y_values = y_values.flatten()                                           #pylint: disable=maybe-no-member
        x_values = x_values.flatten()                                           #pylint: disable=maybe-no-member

        # Do the fit!
        print("Fitting (this can sometimes take a while...)")
        init_resid = self.wave_fit_resid(
            init_mod, orders, x_values, y_values, ydeg=ydeg, xdeg=xdeg)
        bestp = op.leastsq(self.wave_fit_resid, init_mod,
                           args=(orders, x_values, y_values, ydeg, xdeg))
        final_resid = self.wave_fit_resid(
            bestp[0], orders, x_values, y_values, ydeg=ydeg, xdeg=xdeg)
        params = bestp[0].reshape((ydeg + 1, xdeg + 1))
        print(init_resid, final_resid)

        return params

    def spectral_format_with_matrix(self, xmod, wavemod, spatmod=None,
                                    specmod=None, rotmod=None):
        """Create a spectral format, including a detector to slit matrix at
           every point.

           The input parameters are required for this to work and represent
           the polynomial coefficients for second order descriptions of
           how the spectral and spatial scales vary as a function of order
           for each mode as well as a slit rotation indicator.

        The functional form is equivalent to all other models in the spectral
        format. For the 3 new parameters the simplest interpretation of the
        model files is as follows:

        :math:`q_{00}`:  spatial/spectral/rotation scale at the reference order
        :math:`q_{01}`:  variation as a function of orders divided by the number
        of orders
        ... with everything else approximately zero.

        Parameters
        ----------

        xmod: float array
            pixel position model parameters. Used in the spectral format
            function. See documentation there for more details
        wavemod: float array
            pixel position model parameters. Used in the spectral format
            function. See documentation there for more details

        spatmod: (optional) float array
            Parameters from the spatial scale second order polynomial
            describing how the slit image varies in the spatial direction
            as a function of order on the CCD
        specmod: (optional) float array
            Parameters from the spectral scale second order polynomial
            describing how the slit image varies in the spectral direction
            as a function of order on the CCD
        rotmod: (optional) float array
            Parameters from the extra rotation second order polynomial
            describing how the slit image rotation varies
            as a function of order on the CCD

        Returns
        -------
        x: (norders, ny) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel
            and each order (m).
        w: (norders, ny) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (norders, ny) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        matrices: (norders, ny, 2, 2) float array
            2x2 slit rotation matrices, mapping output co-ordinates back
            to the slit.
        """

        if (xmod is None) and (wavemod is None):
            return 'Must provide at least one of xparams or wparams'

        if (spatmod is None) and (specmod is None) and (rotmod is None):
            return 'Must provide at least one of spatmod, specmod or rotmod,\
            otherwise there is no point in running this function.'

        #Get the basic spectral format
        xbase, waves, blaze = self.spectral_format(xparams=xmod,
                                                   wparams=wavemod)
        matrices = np.zeros((xbase.shape[0], xbase.shape[1], 2, 2))
        # Initialise key variables in case models are not supplied.
        amat = np.zeros((2, 2))
        slit_microns_per_det_pix_x = 1.
        slit_microns_per_det_pix_y = 1.
        rotation = 0.

        #Loop through orders
        for order in range(self.m_min, self.m_max + 1):
            mprime = self.m_ref / order - 1
            # Now work out the rotation as a function of pixel along the order
            for yvalue in range(self.szy):
                # Create a matrix where we map input angles to output
                # coordinates.
                # Now obtain the spatial and spectral scales for this order

                # Start with the spatial scale
                if spatmod is not None:
                    # Find the polynomial coefficients for each order.
                    slit_microns_per_det_pix_x = self.evaluate_poly(spatmod,
                                                                    mprime,
                                                                    yvalue)
                if specmod is not None:
                    # Find the polynomial coefficients for each order.
                    slit_microns_per_det_pix_y = self.evaluate_poly(specmod,
                                                                    mprime,
                                                                    yvalue)
                # Then populate the matrix for this location.
                amat[0, 0] = 1.0 / slit_microns_per_det_pix_x
                amat[0, 1] = 0
                amat[1, 0] = 0
                amat[1, 1] = 1.0 / slit_microns_per_det_pix_y

                # Apply an additional rotation matrix. If the simulation was
                # complete, this wouldn't be required.
                # extra_rot should be in radians and come from the third
                # polynomial on the model file.
                if rotmod is not None:
                    rotation = self.evaluate_poly(rotmod, mprime, yvalue)
                r_rad = np.radians(rotation)
                dy_frac = 1. / (xbase.shape[1] / 2.0)
                extra_rot_mat = np.array([[np.cos(r_rad * dy_frac),
                                           np.sin(r_rad * dy_frac)],
                                          [-np.sin(r_rad * dy_frac),
                                           np.cos(r_rad * dy_frac)]])
                amat = np.dot(extra_rot_mat, amat)
                # We actually want the inverse of this (mapping output
                # coordinates back onto the slit.
                matrices[order-self.m_min, yvalue, :, :] = np.linalg.inv(amat)
        return xbase, waves, blaze, matrices

    def spectral_format_with_matrix_fast(self, xmod, wavemod, spatmod=None,
                                    specmod=None, rotmod=None):
        """Create a spectral format, including a detector to slit matrix at
           every point. This is the faster version.

           The input parameters are required for this to work and represent
           the polynomial coefficients for second order descriptions of
           how the spectral and spatial scales vary as a function of order
           for each mode as well as a slit rotation indicator.

        The functional form is equivalent to all other models in the spectral
        format. For the 3 new parameters the simplest interpretation of the
        model files is as follows:

        :math:`q_{00}`:  spatial/spectral/rotation scale at the reference order
        :math:`q_{01}`:  variation as a function of orders divided by the number
        of orders
        ... with everything else approximately zero.

        The explanation of the matrix operations goes as follows:

        We need a 2x2 matrix for every pixel long every order that describes the
        rotation of the slit.
        In order to obtain the final matrix we require a 2x2 diagonal matrix
        containing the inverse of the slit microns per detector pixel scale in
        the spacial and spectral directions as the (0,0) and (1,1) values.
        We then require another 2x2 matrix with the rotation angle built in
        the form:

        \Biggl \lbracket
        { cos(\theta) , sin(\theta) \atop 
          -sin(theta), cos(theta) }
        \rbracket

        These two matrices are then multiplied together with a dot product and
        the desired result is the inverse matrix of the product. 

        However, due to the diagonality of the scale matrix and the time it
        takes to loop through the inversions, we have then explicitly done the
        algebra in this way:

        :math:`( R(\theta) . D )^-1 = D^-1 . R^-1 = D^-1 . R(-\theta)`

        if R is the rotation matrix and D is the scale diagonal matrix

        Parameters
        ----------

        xmod: float array
            pixel position model parameters. Used in the spectral format
            function. See documentation there for more details
        wavemod: float array
            pixel position model parameters. Used in the spectral format
            function. See documentation there for more details

        spatmod: (optional) float array
            Parameters from the spatial scale second order polynomial
            describing how the slit image varies in the spatial direction
            as a function of order on the CCD
        specmod: (optional) float array
            Parameters from the spectral scale second order polynomial
            describing how the slit image varies in the spectral direction
            as a function of order on the CCD
        rotmod: (optional) float array
            Parameters from the extra rotation second order polynomial
            describing how the slit image rotation varies
            as a function of order on the CCD

        Returns
        -------
        x: (norders, ny) float array
            The x-direction pixel co-ordinate corresponding to each y-pixel
            and each order (m).
        w: (norders, ny) float array
            The wavelength co-ordinate corresponding to each y-pixel and each
            order (m).
        blaze: (norders, ny) float array
            The blaze function (pixel flux divided by order center flux)
            corresponding to each y-pixel and each order (m).
        matrices: (norders, ny, 2, 2) float array
            2x2 slit rotation matrices, mapping output co-ordinates back
            to the slit.
        """

        if (xmod is None) and (wavemod is None):
            return 'Must provide at least one of xparams or wparams'

        if (spatmod is None) and (specmod is None) and (rotmod is None):
            return 'Must provide at least one of spatmod, specmod or rotmod,\
            otherwise there is no point in running this function.'

        #Get the basic spectral format
        xbase, waves, blaze = self.spectral_format(xparams=xmod,
                                                   wparams=wavemod)
        matrices = np.zeros((xbase.shape[0], xbase.shape[1], 2, 2))
        # Initialise key variables in case models are not supplied.
        amat = np.zeros((self.szy,2, 2))
        slit_microns_per_det_pix_x = np.ones(self.szy)
        slit_microns_per_det_pix_y = np.ones(self.szy)
        rotation = np.zeros(self.szy)
        yvalue=np.arange(self.szy)
        
        #Loop through orders
        for order in range(self.m_min, self.m_max + 1):
            mprime = self.m_ref / order - 1
            # Now work out the rotation as a function of pixel along the order
            # Create a matrix where we map input angles to output
            # coordinates.
            # Now obtain the spatial and spectral scales for this order

            # Start with the spatial scale
            if spatmod is not None:
                # Find the polynomial coefficients for each order.
                slit_microns_per_det_pix_x = self.evaluate_poly(spatmod,
                                                                mprime,
                                                                yvalue)
            if specmod is not None:
                # Find the polynomial coefficients for each order.
                slit_microns_per_det_pix_y = self.evaluate_poly(specmod,
                                                                mprime,
                                                                yvalue)

            # Apply an additional rotation matrix. If the simulation was
            # complete, this wouldn't be required.
            # extra_rot should be in radians and come from the third
            # polynomial on the model file.
            if rotmod is not None:
                rotation = self.evaluate_poly(rotmod, mprime, yvalue)
            r_rad = np.radians(rotation)
            dy_frac = 1. / (xbase.shape[1] / 2.0)
            extra_rot_mat = np.zeros((self.szy,2, 2))
            # This needs to be done separately because of the
            # matrix structure
            extra_rot_mat[:,0,0] = np.cos(r_rad * dy_frac) * \
                                   slit_microns_per_det_pix_x
            extra_rot_mat[:,0,1] = -np.sin(r_rad * dy_frac) *\
                                   slit_microns_per_det_pix_x
            extra_rot_mat[:,1,0] = np.sin(r_rad * dy_frac) *\
                                   slit_microns_per_det_pix_y
            extra_rot_mat[:,1,1] = np.cos(r_rad * dy_frac) *\
                                   slit_microns_per_det_pix_y
            matrices[order-self.m_min, :, :, :] = extra_rot_mat
        return xbase, waves, blaze, matrices

    def slit_flat_convolve(self, flat, slit_profile=None):
        """Dummy function that would takes a flat field image and a slit profile and
           convolves the two in 2D. Returns result of convolution, which
           should be used for tramline fitting.
        """
        return flat

    def adjust_model(self, data, wparams=None,
                     xparams=None, convolve=True,
                     percentage_variation=10, vary_wrt_max=True):
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
            Vary all parameters intelligently with a scaling of the maximum
            variation, rather than just a percentage of each.

        Returns
        -------
        xparams: 2D float array
             New adjusted x parameters

        """
        # Must provide at least one of the model parameters
        if (xparams is None) and (wparams is None):
            return 'Must provide at least one of xparams or wparams'

        nxbase = data.shape[0]
        # Start by setting up the graphical part
        fig, axx = plt.subplots()
        axcolor = 'lightgoldenrodyellow'
        # Grab the model to be plotted
        x_int, wave_int, blaze_int = self.spectral_format(wparams=wparams,
                                                          xparams=xparams)
        ygrid = np.meshgrid(np.arange(data.shape[1]),
                            np.arange(x_int.shape[0]))[0]
        # The data must be flattened for the sliders to work.
        # Then plot it!
        lplot, = plt.plot(ygrid.flatten()[::10], x_int.flatten()[::10] + nxbase
                          // 2, color='green', linestyle='None', marker='.')

        if convolve:
            data = self.slit_flat_convolve(flat=data)
        # Now over plot the image.
        axx.imshow((data - np.median(data)) / 1e2)

        # Create a second window for sliders.
        slide_fig = plt.figure()

        # This function is executed on each slider change.
        # spectral_format is updated.
        def update(val):
            """ Function used to trigger updates on sliders """
            for i in range(npolys):
                for j in range(polyorder):
                    xparams[i, j] = sliders[i][j].val
            xbase, wave, blaze = self.spectral_format(wparams=wparams,
                                                      xparams=xparams)
            lplot.set_ydata(xbase.flatten()[::10] + nxbase // 2)
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
                frac_xparams[i] = np.max(
                    frac_xparams[-1]) / (nxbase / 2.0)**(npolys - 1 - i)
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
            """Function for the button tasks"""
            plt.close('all')
            return xparams

        button.on_clicked(submit)

        plt.show()
        return xparams


    def adjust_model_with_matrices(self, data, xparams,
                                   wparams, spatmod, specmod,
                                   rotmod, percentage_variation=10,
                                   vary_wrt_max=True):
        """Function that uses matplotlib slider widgets to adjust a polynomial
        model overlaid on top of an image. In practice this will be
        overlaid on top of the result of convolving the flat with a slit
        profile in 2D for the xmod only, and then an arc frame for all other 
        parameters. User MUST provide an initial model for all 5 aspects.

        This function actually works by displaying the result of the model on
        top of an image of some sort, while bringing up other windows with 
        either sliders or value boxes for the polynomial parameters and include
        update buttons for the user to play around with values. This should 
        be used by the engineering team to manually work out an initial model
        before running the fitting function. 

        Parameters
        ----------
        data: float array
            an array containing data to be used as a visual comparison of the
            model
        xparams: float array
            2D array containing the initial wavelength
            model parameters.
        wparams: float array
            2D array containing the initial order location model parameters.
        spatmod: float array
            2D array containing the initial spatial direction
            slit model parameters.
        specmod: float array
            2D array containing the initial spectral direction
            slit model.
        rotmod: float array
            2D array containing the initial slit rotation model parameters.
        percentage_variation: int
            How much should the percentage adjustment in each bin as a function
            of the parameter. Default is 10%
        vary_wrt_max: bool
            Vary all parameters intelligently with a scaling of the maximum
            variation, rather than just a percentage of each.

        Returns
        -------
        Either a success or an error message.

        """
        

        nxbase = data.shape[0]
        # Start by setting up the graphical part
        fig, axx = plt.subplots()
        axcolor = 'lightgoldenrodyellow'
        # Grab the model to be plotted
        x_int, wave_int, blaze_int = self.spectral_format(wparams=wparams,
                                                          xparams=xparams)
        ygrid = np.meshgrid(np.arange(data.shape[1]),
                            np.arange(x_int.shape[0]))[0]
        # The data must be flattened for the sliders to work.
        # Then plot it!
        lplot, = plt.plot(ygrid.flatten()[::10], x_int.flatten()[::10] + nxbase
                          // 2, color='green', linestyle='None', marker='.')

        if convolve:
            data = self.slit_flat_convolve(flat=data)
        # Now over plot the image.
        axx.imshow((data - np.median(data)) / 1e2)

        # Create a second window for sliders.
        slide_fig = plt.figure()

        # This function is executed on each slider change.
        # spectral_format is updated.
        def update(val):
            """ Function used to trigger updates on sliders """
            for i in range(npolys):
                for j in range(polyorder):
                    xparams[i, j] = sliders[i][j].val
            xbase, wave, blaze = self.spectral_format(wparams=wparams,
                                                      xparams=xparams)
            lplot.set_ydata(xbase.flatten()[::10] + nxbase // 2)
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
                frac_xparams[i] = np.max(
                    frac_xparams[-1]) / (nxbase / 2.0)**(npolys - 1 - i)
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
            """Function for the button tasks"""
            plt.close('all')
            return xparams

        button.on_clicked(submit)

        plt.show()
        return xparams
