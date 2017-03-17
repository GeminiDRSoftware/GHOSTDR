
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

        if not slit_profile:
            # At this point create a slit profile
            # Create a x baseline for convolution and assume a FWHM for profile
            xbase = flat.shape[0]
            profilex = np.arange(xbase) - xbase // 2
            # Now create a model of the slit profile
            mod_slit = np.zeros(xbase)
            if self.mode == 'high':
                nfibers = 26
            else:
                nfibers = self.nlenslets

            for i in range(-(nfibers // 2), -(nfibers // 2) + nfibers):
                mod_slit += np.exp(-(profilex - i * fiber_separation)**2 /
                                   2.0 / profile_sigma**2)
        else:
            mod_slit = slit_profile

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
