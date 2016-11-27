""" This is a simple python script to test the simulator code. """

import math
import numpy as np
import pyghost
import pylab as plt

def run(nbias=3, ndark=3, nflat=3, crplane=False):
    """ The function that runs the test. """

    # Create the two arms
    blue = pyghost.Arm('blue')
    red = pyghost.Arm('red')

    # Create a blank spectrum (used for the bias, dark, and sky)
    blank = np.array([[0.1, 1.0], [0.0, 0.0]])

    # Create a perfectly flat spectrum (used for the flat); need to replace
    # this with a better GCAL spectrum
    flat = np.array([[0.1, 1.0], [100.0, 100.0]])

    # Create a ThAr spectrum (used for the arc)
    thar = pyghost.thar_spectrum()

    # We read out through 4 amps
    namps = [2, 2]

    # Our default bias level
    bias_level = 100

    # Our default overscan
    overscan = 32

    # Our default duration
    duration = 100

    noise = {}
    scaling = {}

    for arm in (blue, red):
        # Generate some noise with regular 50 and 60 Hz frequencies.
        fnoise = pyghost.split_image(arm.simulate_frequency_noise([50, 60], 0, 5), namps, return_headers=False)

        # Generate a gradient over the whole detector
        gnoise = pyghost.split_image(arm.simulate_gradient(math.pi/4, 0, 10), namps, return_headers=False)

        # The additive noise is the sum of these two
        noise[arm.arm] = fnoise + gnoise

        # Generate the pixel-to-pixel variation image, with 1% rms error
        scaling[arm.arm] = pyghost.split_image(arm.simulate_flatfield(1.0, 0.01), namps, return_headers=False)

    # This object captures the details of the detectors
    ghost = pyghost.Ghost(rnoise=3.0, gain=1.0, namps=namps, overscan=overscan, bias_level=bias_level,
                          additive_noise=noise, scaling=scaling)

    ## This produces a bias with the above noise
    for i in range(1, nbias+1):
        ghost.simulate_observation(duration=0.0, output_prefix='bias_{0:d}_'.format(i),
                                   spectrum_in=blank, use_thar=False,
                                   add_sky=False, obstype='BIAS',
                                   write_crplane=crplane, data_label=i)

    ## This produces a dark frame
    for i in range(1, ndark+1):
        ghost.simulate_observation(duration=duration,
                                   output_prefix='dark'+str(duration)+'_{0:d}_'.format(i),
                                   use_thar=False, spectrum_in=blank, add_sky=False,
                                   obstype='DARK', write_crplane=crplane, data_label=i)

    for mode in ('std', 'high'):
        # This (should) produce a GCAL flat frame
        for i in range(1, nflat+1):
            ghost.simulate_observation(duration=duration,
                                       output_prefix='flat'+str(duration)+'_'+mode+'_{0:d}_'.format(i),
                                       use_thar=False, spectrum_in=flat,
                                       add_sky=False, mode=mode,
                                       flatlamp=True, obstype='FLAT',
                                       write_crplane=crplane, data_label=i)

        # This produces an arc frame
        ghost.simulate_observation(duration=duration,
                                   output_prefix='arc'+str(duration)+'_'+mode+'_',
                                   use_thar=False, spectrum_in=thar, add_sky=False,
                                   mode=mode, flatlamp=True, obstype='ARC',
                                   write_crplane=crplane)

        # This produces a sky frame
        ghost.simulate_observation(duration=duration,
                                   output_prefix='sky'+str(duration)+'_'+mode+'_',
                                   use_thar=False, spectrum_in=blank,
                                   add_sky=True, mode=mode, obstype='SKY',
                                   write_crplane=crplane)

        # Make the slit-viewing flux a bit random, to simulate clouds
        svfp = np.random.randn(100)
        # Shift and normalise to [0, 1]
        svfp -= svfp.min()
        svfp /= svfp.max()
        sv_flux_profile = (np.arange(len(svfp)), svfp)

        # This produces an object frame, using the default object spectrum
        # and 0.5 arcsec seeing
        ghost.simulate_observation(duration=duration,
                                   output_prefix='obj'+str(duration)+'_0.5_'+mode+'_',
                                   use_thar=True, add_sky=True, mode=mode,
                                   obstype='OBJECT', seeing=0.5,
                                   write_crplane=crplane, sv_flux_profile=sv_flux_profile)

        # This produces an object frame, using the default object spectrum
        # and 1.0 arcsec seeing
        ghost.simulate_observation(duration=duration,
                                   output_prefix='obj'+str(duration)+'_1.0_'+mode+'_',
                                   use_thar=True, add_sky=True, mode=mode,
                                   obstype='OBJECT', seeing=1.0,
                                   write_crplane=crplane, sv_flux_profile=sv_flux_profile)

if __name__ == "__main__":
    run()
