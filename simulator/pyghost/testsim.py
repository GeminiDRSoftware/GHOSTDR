""" This is a simple python script to test the simulator code. """

import numpy as np
import pyghost


def run_test():
    """ The function that runs the test. """

    # Create the two arms
    blue = pyghost.Arm('blue')
    red = pyghost.Arm('red')

    # Create a blank spectrum (used for the bias, dark, and sky)
    blank = np.array([[0.1, 1.0], [0.0, 0.0]])

    # Create a ThAr spectrum (used for the arc)
    thar = blue.thar_spectrum()

    for arm in (red, blue):
        # This should produce a completely blank frame
        arm.simulate_frame(duration=0.0, output_prefix='zero_', use_thar=False,
                           rnoise=0.0, spectrum=blank, gain=1, namps=[2, 2],
                           overscan=32, add_sky=False, bias_level=0)

        # This produces a realistic looking bias
        arm.simulate_frame(duration=0.0, output_prefix='bias_', use_thar=False,
                           rnoise=3.0, spectrum=blank, gain=1, namps=[2, 2],
                           overscan=32, add_sky=False, bias_level=100)

        # This produces an arc frame
        # FIXME we need an option to say the input spectrum is
        # evenly illuminating the fibers
        arm.simulate_frame(duration=300.0, output_prefix='arc300_',
                           use_thar=False, rnoise=3.0, spectrum=thar, gain=1,
                           namps=[2, 2], overscan=32, add_sky=False,
                           bias_level=100)

        # This produces a dark frame
        arm.simulate_frame(duration=300.0, output_prefix='dark300_',
                           use_thar=False, rnoise=3.0, spectrum=blank, gain=1,
                           namps=[2, 2], overscan=32, add_sky=False,
                           bias_level=100)

        # This produces a sky frame
        arm.simulate_frame(duration=300.0, output_prefix='sky300_',
                           use_thar=False, rnoise=3.0, spectrum=blank, gain=1,
                           namps=[2, 2], overscan=32, add_sky=True,
                           bias_level=100)

        # This produces an object frame, using the default object spectrum
        arm.simulate_frame(duration=300.0, output_prefix='obj300_',
                           use_thar=True, rnoise=3.0, gain=1, namps=[2, 2],
                           overscan=32, add_sky=True, bias_level=100)


if __name__ == "__main__":
    run_test()
