""" This is a simple python script to test the simulator code. """

import math
import numpy as np
import pyghost

def run(nbias=3,ndark=3,nflat=3,crplane=False):
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

    for arm in (red, blue):
        # This should produce a completely blank frame
        #arm.simulate_frame(duration=0.0, output_prefix='zero_', use_thar=False,
        #                   rnoise=0.0, spectrum=blank, gain=1, namps=namps,
        #                   overscan=0, add_sky=False, bias_level=0,
        #                   obstype='BIAS')

        # This produces a realistic looking bias, without pattern noise
        #arm.simulate_frame(duration=0.0, output_prefix='bias_', use_thar=False,
        #                   rnoise=3.0, spectrum=blank, gain=1, namps=namps,
        #                   overscan=overscan, add_sky=False,
        #                   bias_level=bias_level, obstype='BIAS')

        # Generate some noise with regular 50 and 60 Hz frequencies.
        fnoise = pyghost.split_image(arm.simulate_frequency_noise([50, 60], 0, 5), namps, return_headers=False)

        # Generate a gradient over the whole detector
        gnoise = pyghost.split_image(arm.simulate_gradient(math.pi/4, 0, 10), namps, return_headers=False)

        # The additive noise is the sum of these two
        noise = fnoise + gnoise

        # Generate the pixel-to-pixel variation image, with 1% rms error
        scaling = pyghost.split_image(arm.simulate_flatfield(1.0, 0.01), namps, return_headers=False)

        # This produces a bias with the above noise
        for i in range(nbias):
            arm.simulate_frame(duration=0.0, output_prefix='bias_{0:d}_'.format(i), use_thar=False,
                               rnoise=3.0, spectrum_in=blank, gain=1, namps=namps,
                               overscan=overscan, add_sky=False, bias_level=bias_level,
                               obstype='BIAS', additive_noise=noise, scaling=scaling,
                               write_crplane=crplane)

        # This produces a dark frame
        for i in range(ndark):
            arm.simulate_frame(duration=duration,
                               output_prefix='dark'+str(duration)+'_{0:d}_'.format(i),
                               use_thar=False, rnoise=3.0, spectrum_in=blank, gain=1,
                               namps=namps, overscan=overscan, add_sky=False,
                               bias_level=bias_level, obstype='DARK',
                               additive_noise=noise, scaling=scaling,
                               write_crplane=crplane)

        for mode in ('std', 'high'):
            # This (should) produce a GCAL flat frame
            for i in range(nflat):
                arm.simulate_frame(duration=duration,
                                   output_prefix='flat'+str(duration)+'_'+mode+'_{0:d}_'.format(i),
                                   use_thar=False, rnoise=0.0, spectrum_in=flat, gain=1,
                                   namps=namps, overscan=overscan, add_sky=False, mode=mode,
                                   bias_level=bias_level, flatlamp=True, obstype='FLAT',
                                   additive_noise=noise, scaling=scaling,
                                   write_crplane=crplane)

            # This produces an arc frame
            arm.simulate_frame(duration=duration,
                               output_prefix='arc'+str(duration)+'_'+mode+'_',
                               use_thar=False, rnoise=3.0, spectrum_in=thar, gain=1,
                               namps=namps, overscan=overscan, add_sky=False,
                               mode=mode, bias_level=bias_level, flatlamp=True,
                               obstype='ARC', additive_noise=noise, scaling=scaling,
                               write_crplane=crplane)

            # This produces a sky frame
            arm.simulate_frame(duration=duration,
                               output_prefix='sky'+str(duration)+'_'+mode+'_',
                               use_thar=False, rnoise=3.0, spectrum_in=blank, gain=1,
                               namps=namps, overscan=overscan, add_sky=True,
                               mode=mode, bias_level=bias_level, obstype='SKY',
                               additive_noise=noise, scaling=scaling,
                               write_crplane=crplane)

            # This produces an object frame, using the default object spectrum
            # and 0.5 arcsec seeing
            arm.simulate_frame(duration=duration,
                               output_prefix='obj'+str(duration)+'_0.5_'+mode+'_',
                               use_thar=True, rnoise=3.0, gain=1, namps=namps,
                               overscan=overscan, add_sky=True, mode=mode,
                               bias_level=bias_level, obstype='OBJECT',
                               additive_noise=noise, scaling=scaling, seeing=0.5,
                               write_crplane=crplane)

            # This produces an object frame, using the default object spectrum
            # and 1.0 arcsec seeing
            arm.simulate_frame(duration=duration,
                               output_prefix='obj'+str(duration)+'_1.0_'+mode+'_',
                               use_thar=True, rnoise=3.0, gain=1, namps=namps,
                               overscan=overscan, add_sky=True, mode=mode,
                               bias_level=bias_level, obstype='OBJECT',
                               additive_noise=noise, scaling=scaling, seeing=1.0,
                               write_crplane=crplane)

    #Now combine the slit images.
    for rfn, bfn in zip(red.slit_fns, blue.slit_fns):
        pyghost.ghostsim.combine_slitviewer_files([rfn,bfn])

if __name__ == "__main__":
    run()
