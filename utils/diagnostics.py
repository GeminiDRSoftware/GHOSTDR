"""A script to do diagnostic checks on reduced ghost data.

It displays several windows showing how good the fits are and a few other things

This script should be ran from within the finished reduction folder, which 
should contain the 'calibrations/' directory that are needed and finished 
reduction files (extracted profiles/barycentric corrected)
"""

from __future__ import division, print_function
import numpy as np
from ghostdr import polyfit
import glob, os, sys
import astropy.io.fits as pyfits
import ghostdr.ghost.lookups as lookups

# Start by finding where the lookups are and save that location as a variable
lookups_path = os.path.dirname( os.path.abspath(lookups.__file__) )


def thar_spectrum(linefile):
    """Calculates a ThAr spectrum. Note that the flux scaling here is roughly
    correct for the lamp with no neutral density. From the simulator.

    Parameters
    ----------

    linefile: string
        Path to the arcline file


    Returns
    -------
    wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
    """

    thar = np.loadtxt(linefile, usecols=[0, 1, 2])
    # Create a fixed wavelength scale evenly spaced in log.
    thar_wave = 3600 * np.exp(np.arange(5e5) / 5e5)
    thar_flux = np.zeros(int(5e5))
    # NB This is *not* perfect: we just find the nearest index of the
    # wavelength scale corresponding to each Th/Ar line.
    wave_ix = (np.log(thar[:, 1] / 3600) * 5e5).astype(int)                     
    wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5 - 1).astype(int)           
    thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))                        
    thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                            mode='same')

    thar_flux /= np.max(thar_flux) / 3e5
    return np.array([thar_wave, thar_flux])



# Let's start by checking the fits. We use the same method as the slider adjust
# and let the user check things individually.


flat_list = glob.glob('calibrations/processed_flat/*flat.fits')
arc_list = glob.glob('calibrations/processed_arc/*arc.fits')

# Now cycle through available modes.

for mode in ['high', 'std']:
    for cam in ['blue', 'red']:
        ghost = polyfit.ghost.GhostArm(cam, mode=mode)
        
