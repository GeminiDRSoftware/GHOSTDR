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
from __future__ import division, print_function
import numpy as np
from astrodata_GHOST.polyfit.polyspect import Polyspect

class GhostArm(Polyspect):
    """A class for each arm of the spectrograph. The initialisation
    function takes a series of strings representing the configuration.
    It can be "red" or "blue" for the arm (first string),
    and "std" or "high" for the mode (second string).
    """

    def __init__(self, arm='blue', mode='std'):
        """Initialisation function that sets all the mode specific parameters
        related to each configuration of the spectrograph.

        Attributes
        ----------
        
        arm: str
            Which arm of the GHOST spectrograph is to be initialized.
        spect: str
            Which spectrograph in usage.
        lenslet_high_size: int
            Lenslet flat-to-flat in microns for high mode
        lenslet_std_size: int
            Lenslet flat-to-flat in microns for standard mode
        mode: str
            Resolution mode.
        nlenslets: int
            number of lenslets of the IFU
        """
        if arm == 'red':
            Polyspect.__init__(self, m_ref=50, szx=6144, szy=6160, m_min=34,
                               m_max=67, transpose=True)
        elif arm == 'blue':
            Polyspect.__init__(self, m_ref=80, szx=4096, szy=4112, m_min=63,
                               m_max=95, transpose=True)
        else:
            print("Unknown spectrograph arm!")
            raise UserWarning
        # A lot of these parameters are yet unused.
        # These are a legacy of the original simulator and are left here
        # because they may become useful in the future.
        self.spect = 'ghost'
        self.arm = arm
        self.lenslet_high_size = 118.0  # Lenslet flat-to-flat in microns
        self.lenslet_std_size = 197.0  # Lenslet flat-to-flat in microns
        self.mode = mode

        # Now we determine the number of fibers based on mode.
        if mode == 'high':
            self.nlenslets = 28
        elif mode == 'std':
            self.nlenslets = 17
        else:
            print("Unknown mode!")
            raise UserWarning

    def bin_data(self,data, binning=[2,2]):
        """ Generic Function used to create a binned equivalent of a spectrograph
        image array for the purposes of equivalent extraction. 

        Parameters
        ----------
        data: :obj:`numpy.ndarray'
            The (unbinned) data to be binned
        binning: list
            A two element list with the binning factors

        Returns
        -------
        binned_array: :obj:`numpy.ndarray'
            The binned data
        """
        if data.shape != (self.szx,self.szy):
            raise UserWarning('Input data for binning is not in the expected\
            format')
        
        rows = binning[0]
        cols = binning[1]
        binned_array = data.reshape(int(data.shape[0]/rows),
                                    rows,
                                    int(data.shape[1]/cols),
                                    cols).\
                                    sum(axis=1).sum(axis=2)
        return binned_array

        
    def slit_flat_convolve(self, flat, slit_profile=None, spatpars=None,\
        microns_pix=None, xpars=None, num_conv=3):
        """Function that takes a flat field image and a slit profile and
           convolves the two in 2D. Returns result of convolution, which
           should be used for tramline fitting.

        Parameters
        ----------
        flat: :obj:`numpy.ndarray`
            A flat field image from the spectrograph
            
        slit_profile: :obj:`numpy.ndarray`, optional
            A slit profile as a 1D array with the slit profile fiber amplitudes.
            If none is supplied this function will assume identical fibers and
            create one to be used in the convolution based on default parameters
            specified in the ghost class.
            
        spatpars: :obj:`numpy.ndarray`, optional
            The 2D polynomial parameters for the slit spatial scale. 
            Required if slit_profile is not None.
        
        microns_pix: float, optional
            The slit scale in microns per pixel
            Required if slit_profile is not None.
            
        xpars:  :obj:`numpy.ndarray`, optional
            The 2D polynomial parameters for the x (along-slit) coordinate.
            Required if slit_profile is not None.
            
        num_conv: int, optional, optional
            The number of different convolution functions to use for different
            orders.
            The final convolved profile is an interpolation between these.

        Returns
        -------
        flat_conv: :obj:`numpy.ndarray`
            The convolved 2D array.
        """
        #FIXME: Error checking of inputs is needed here
        #TODO: Based on test of speed, the convolution code with an input slit_profile
        #      could go order-by-order.

        if self.arm == 'red':
            # Now put in the default fiber profile parameters for each mode.
            # These are used by the convolution function on polyspect
            # These were determined based on visual correspondence with
            # simulated data and may need to be revised once we have real
            # data. The same applies to the blue arm parameters.
            if self.mode == 'std':
                fiber_separation = 4.15
                profile_sigma = 1.1
            elif self.mode == 'high':
                fiber_separation = 2.49
                profile_sigma = 0.7
        elif self.arm == 'blue':
            # Additional slit rotation accross an order needed to match Zemax.
            #self.extra_rot = 2.0

            # Now put in the default fiber profile parameters for each mode.
            # These are used by the convolution function on polyspect
            if self.mode == 'std':
                fiber_separation = 3.97
                profile_sigma = 1.1
            elif self.mode == 'high':
                fiber_separation = 2.53
                profile_sigma = 0.7
        else:
            print("Unknown spectrograph arm!")
            raise UserWarning

        # Fourier transform the flat for convolution
        im_fft = np.fft.rfft(flat, axis=0)
        
        # Create a x baseline for convolution 
        xbase = flat.shape[0]
        profilex = np.arange(xbase) - xbase // 2

        #This is the original code which is based on the fixed fiber_separation
        #defined above. 
        if slit_profile is None:
            flat_conv = np.zeros_like(im_fft)
            
            # At this point create a slit profile
            
            # Now create a model of the slit profile
            mod_slit = np.zeros(xbase)
            if self.mode == 'high':
                nfibers = 26
            else:
                nfibers = self.nlenslets

            for i in range(-(nfibers // 2), -(nfibers // 2) + nfibers):
                mod_slit += np.exp(-(profilex - i * fiber_separation)**2 /
                                   2.0 / profile_sigma**2)
            # Normalise the slit model and fourier transform for convolution
            mod_slit /= np.sum(mod_slit)
            mod_slit_ft = np.fft.rfft(np.fft.fftshift(mod_slit))
            
            # Now convolved in 2D
            for i in range(im_fft.shape[1]):
                flat_conv[:, i] = im_fft[:, i] * mod_slit_ft
            
            #Now inverse transform.
            flat_conv = np.fft.irfft(flat_conv, axis=0)
            
        else:
            flat_conv = np.zeros_like(flat)
            flat_conv_cube = np.zeros( (num_conv, flat.shape[0], flat.shape[1]) )
            
            #Our orders that we'll evaluate the spatial scale at:
            orders = np.linspace(self.m_min, self.m_max, num_conv).astype(int)
            mprimes = self.m_ref / orders - 1
            y_values = np.arange(self.szy)
            
            #The slit coordinate in microns
            slit_coord = (np.arange(len(slit_profile)) -
                          len(slit_profile)//2) * microns_pix

            x_map = np.empty( (len(mprimes), self.szy) )
            
            # Now convolved in 2D
            for j, mprime in enumerate(mprimes):
            #The spatial scales
                spat_scale = self.evaluate_poly(spatpars)[orders[j]-self.m_min]
                
                #The x pixel values, just for this order
                x_map[j] = self.evaluate_poly(xpars)[orders[j]-self.m_min]
                for i in range(im_fft.shape[1]):
                    #Create the slit model.
                    mod_slit = np.interp(profilex*spat_scale[i], slit_coord, slit_profile)
                    
                    # Normalise the slit model and Fourier transform for convolution
                    mod_slit /= np.sum(mod_slit)
                    mod_slit_ft = np.fft.rfft(np.fft.fftshift(mod_slit))
                    #FIXME: Remove num_conv on next line and see if it makes a difference!
                    flat_conv_cube[j, :, i] = np.fft.irfft((im_fft[:, i] * mod_slit_ft)/num_conv)
            
            #Work through every y coordinate and interpolate between the convolutions
            #with the different slit profiles.
            x_ix = np.arange(flat.shape[0])- flat.shape[0]//2

            #Create an m index, and reverse x_map if needed.
            #FIXME: This assumes a minimum size of x_map which should be checked above,
            #i.e. mprimes has 2 or more elements.
            m_map_ix = np.arange(len(mprimes))
            if x_map[1,0] < x_map[0,0]:
                m_map_ix = m_map_ix[::-1]
                x_map = x_map[::-1]
            for i in range(im_fft.shape[1]):
                m_ix_for_interp = np.interp(x_ix, x_map[:,i], m_map_ix)
                m_ix_for_interp = np.minimum(m_ix_for_interp, len(mprimes)-1-1e-6)
                m_ix_for_interp = np.maximum(m_ix_for_interp, 0)
                m_ix_lo = np.int16(m_ix_for_interp)
                m_ix_hi = m_ix_lo+1
                m_ix_frac = m_ix_for_interp - m_ix_lo
                for j in range(len(mprimes)):
                    weight = (m_ix_lo==j) * (1-m_ix_frac) + (m_ix_hi==j) * m_ix_frac
                    flat_conv[:,i] += weight * flat_conv_cube[j, :, i]
            
        return flat_conv
