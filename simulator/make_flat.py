# script that can be used to create an image.
import pyghost

g=pyghost.ghostsim.Arm('red',cosmics=False, crplane=False,split=True)
flat=g.get_solar_spectrum()
flat[1]=120.0
g.simulate_frame(duration=120.0, output_prefix='flathigh_', spectrum_in=flat,
                 bias_level=0, overscan=0, namps=[1, 1],
                 use_thar=False, res='high',
                 add_sky=False, return_image=False, thar_flatlamp=False,
                 flatlamp=True, obstype='FLAT', additive_noise=None,
                 scaling=None, seeing=0.8)
