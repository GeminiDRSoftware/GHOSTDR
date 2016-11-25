# script that can be used to create an image.
import pyghost

g=pyghost.ghostsim.Arm('red')
flat=g.get_solar_spectrum()
flat[1]=120.0
g.simulate_frame(duration=120.0, output_prefix='flatstd_', spectrum_in=flat,                     bias_level=0, overscan=0, namps=[1, 1],
                       use_thar=False, mode='std', add_cosmic=False,
                       add_sky=False, return_image=False, thar_flatlamp=False,                    flatlamp=True, obstype=None, additive_noise=None,
                       scaling=None, seeing=0.8, write_crplane=False)
