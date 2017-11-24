""" This is a simple python script to test the simulator code. """

import math
import numpy as np
import pyghost

def run(nbias=3, ndark=3, nflat=3, cosmics=True, crplane=False, hpplane=False,
        split=True, check=False):
    """ The function that runs the test. """

    # Create the two arms
    blue = pyghost.Arm('blue', cosmics=cosmics, crplane=crplane,
                       hpplane=hpplane, split=split, check=check)
    red = pyghost.Arm('red', cosmics=cosmics, crplane=crplane,
                      hpplane=hpplane, split=split, check=check)

    # Create a blank spectrum (used for the bias, dark, and sky)
    blank = np.array([[0.1, 1.2], [0.0, 0.0]])

    # Create a perfectly flat spectrum (used for the flat); need to replace
    # this with a better GCAL spectrum
    flat = np.array([[0.1, 1.2], [100.0, 100.0]])

    # Create a ThAr spectrum (used for the arc)
    thar = pyghost.thar_spectrum()

    # We read out through 4 amps
    namps = [2, 2]

    # Our default bias level
    bias_level = 100

    # Our default overscan
    oscan = 32

    # Our default duration (pick something not a multiple of the default slit
    # viewer duration [10] so that partial science/slitviewer exposure overlaps
    # can be tested)
    duration = 95

    noise = {}
    scaling = {}

    binmodes = [(1, 1), (1, 2), (1, 8), (2, 4), (2, 8)]

    #Standards list
    stds = ['hd160617.fits','hd180609.fits','hd200654.fits']
    
    for arm in (blue, red):
        # Generate the pixel-to-pixel variation image, with 1% rms error
        scale = arm.simulate_flatfield(1.0, 0.01, oscan)
        scaling[arm.arm] = pyghost.split_image(scale, namps)

        noise[arm.arm] = {}
        for bmode in binmodes:
            # Generate some noise with regular 50 and 60 Hz frequencies.
            fnoise = arm.simulate_frequency_noise([50, 60], 0, 5, bmode, oscan)

            # Generate a gradient over the whole detector
            gnoise = arm.simulate_gradient(math.pi/4, 0, 10, bmode, oscan)

            # The additive noise is the sum of these two
            noise[arm.arm][bmode] = pyghost.split_image(fnoise+gnoise, namps)

    # This object captures the details of the detectors
    ghost = pyghost.Ghost(
        rnoise=3.0, gain=1.0, namps=namps, overscan=oscan,
        bias_level=bias_level, additive_noise=noise, scaling=scaling,
        cosmics=cosmics, crplane=crplane, hpplane=hpplane, split=split, check=check)

    target_binmode = (1, 1)  # (1, 2)

    # This produces a bias with the above noise
    for i in range(1, nbias+1):
        ghost.simulate_observation(
            duration=0.0, output_prefix='bias_'+str(i)+'_',
            spectrum_in=blank, use_thar=False, add_sky=False,
            obstype='BIAS', data_label=i, binmode=target_binmode)

    # This produces a dark frame
    for i in range(1, ndark+1):
        ghost.simulate_observation(
            duration=duration, output_prefix='dark'+str(duration)+'_'+str(i)+'_',
            use_thar=False, spectrum_in=blank, add_sky=False,
            obstype='DARK', data_label=i, binmode=target_binmode)

    for res in ('std', 'high'):
        # This (should) produce a GCAL flat frame
        for i in range(1, nflat+1):
            ghost.simulate_observation(
                duration=duration, output_prefix='flat'+str(duration)+'_'+res +
                '_'+str(i)+'_', use_thar=False, spectrum_in=flat,
                add_sky=False, res=res, flatlamp=True, obstype='FLAT',
                data_label=i, binmode=target_binmode)

        # This produces an arc frame
        ghost.simulate_observation(
            duration=duration, output_prefix='arc'+str(duration)+'_'+res+'_',
            use_thar=False, spectrum_in=thar, add_sky=False, res=res,
            flatlamp=True, obstype='ARC', binmode=target_binmode)

        # This produces a sky frame
        ghost.simulate_observation(
            duration=duration, output_prefix='sky'+str(duration)+'_'+res+'_',
            use_thar=False, spectrum_in=blank, add_sky=True, res=res,
             obstype='SKY', binmode=target_binmode)

        # Make the slit-viewing flux a bit random, to simulate clouds
        svfp = np.random.randn(100)
        # Shift and normalise to [0, 1]
        svfp -= svfp.min()
        svfp /= svfp.max()
        sv_flux_profile = (np.arange(len(svfp)), svfp)

        # Make standard star observations 
        for std in stds:
            spectrum = ghost.get_standard_spectrum(std=std)
            ghost.simulate_observation(
                duration=duration, output_prefix='standard'+str(duration)+'_' +
                std[:-5]+'_'+res+'_', use_thar=True, spectrum_in = spectrum,
                add_sky=True, res=res, obstype='STANDARD', objname=std[:-5],
                data_label=1, sv_flux_profile=sv_flux_profile,
                binmode=target_binmode)

        for i, seeing in enumerate((0.5, 1.0), start=1):  # in arcsecs
            ghost.simulate_observation(
                duration=duration, output_prefix='obj'+str(duration)+'_' +
                str(seeing)+'_'+res+'_', use_thar=True, add_sky=True, res=res,
                obstype='OBJECT', seeing=seeing, data_label=i,
                sv_flux_profile=sv_flux_profile, binmode=target_binmode)

if __name__ == "__main__":
    run()
