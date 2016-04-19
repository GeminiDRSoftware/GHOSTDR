""" NB: Both Sampling and resolution are lowest at the low wavelength ends of each
order. """

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from utils import * #!!! This should be an astro-optics import
plt.ion()

infile = 'HR_slit.txt'
scaling_for_tolerances = 85.0/82.6
outfn = 'HR_Slit_convolved_18.txt'

#infile = 'SR_slit.txt'
#scaling_for_tolerances = 142.0/141.0

include_slit_tilt=True

aberration_sigma = 0.3058 #For <870nm
#aberration_sigma = 0.4587 #For >870nm

#0.3058 corresponds to 80% ensquared energy in 1 pixel and 0.4587 corresponds to 80% ensquared energy in 1.5 pixels.
hex_apsize = 197.0/1e3 #This is a very minor perturbation on Ross's profile.
pixsize_slit_units = 142.0/1e3/2.7
res_slitwidth_product = 51000*0.142


#-----

mm_pix = 0.001
sz = 512

x = (np.arange(sz)-sz//2)*mm_pix
xy = np.meshgrid(x,x)
rr = np.sqrt(xy[0]**2 + xy[1]**2)
r_I = np.loadtxt(infile)
r_I[:,0] *= scaling_for_tolerances
II = np.interp(rr, r_I[:,0], r_I[:,1])
II += np.interp(rr, -r_I[:,0], r_I[:,1])
hex_aperture = hexagon(sz, hex_apsize/mm_pix)

prof_in = np.sum(II*hex_aperture,axis=1)
prof_in /= np.max(prof_in)
prof_pix = np.zeros(sz)
prof_pix[sz//2 - int(pixsize_slit_units/mm_pix/2):sz//2 + int(pixsize_slit_units/mm_pix/2)]=1.0
prof_ab = np.exp(- (x/(aberration_sigma*pixsize_slit_units))**2 /2.0)

#Shift to zero
prof_ab = np.roll(prof_ab,sz//2)
prof_pix = np.roll(prof_pix,sz//2)

if include_slit_tilt:
    final_prof = np.fft.irfft(np.fft.rfft(prof_pix)*np.fft.rfft(prof_pix)*np.fft.rfft(prof_in)*np.fft.rfft(prof_ab))
else:
    final_prof = np.fft.irfft(np.fft.rfft(prof_pix)*np.fft.rfft(prof_in)*np.fft.rfft(prof_ab))
final_prof /= np.max(final_prof)

plt.clf()
plt.plot(x, prof_pix)
plt.plot(x,prof_in)
plt.plot(x,prof_ab)
plt.plot(x,final_prof)

fwhm = np.interp(0.5, final_prof[:sz//2:-1], x[:sz//2:-1]) - \
       np.interp(0.5, final_prof[:sz//2], x[:sz//2])
R = res_slitwidth_product/fwhm
print("Resolution: {0:8.1f}".format(R))

if len(outfn)>0:
    np.savetxt(outfn, np.array([x/pixsize_slit_units,final_prof]).T, fmt='%.4f')