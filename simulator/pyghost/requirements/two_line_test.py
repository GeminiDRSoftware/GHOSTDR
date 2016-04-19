from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op

plt.ion()

prof_fn = 'HR_Slit_convolved_18.txt'
km_sec_pix = 2.0 * (2.0/1.8)
fwhm = 1.8 #In pix

#prof_fn = 'HR_Slit_convolved_20.txt'
#km_sec_pix = 2.0 
#fwhm = 2.0 #In pix

line_width = 2.0 #In km/s.
nphase = 1000
line_sep = 6.0
line_rat = 0.5
snr = 100
#-------------------------

def two_line_resid(params,x_subsamp, gauss_prof, x_data, data, sig):
    """This function is designed for scipy.optimize.leastsq
    Parameters
    ----------
    params: float array
        x1, amp1, x2, amp2
    """
    model_two_lines =  params[1]*np.interp(x_data - params[0],x_subsamp,gauss_prof)
    model_two_lines += params[3]*np.interp(x_data - params[2],x_subsamp,gauss_prof)
    return (model_two_lines - data)/sig

#Create a realistic line profile.
prof = np.loadtxt(prof_fn)
x_subsamp = prof[:,0]
g = np.exp(-x_subsamp**2/2.0/(line_width/km_sec_pix)**2)
g /= np.sum(g)
line = np.fft.fftshift(np.fft.irfft(np.fft.rfft(g)*np.fft.rfft(prof[:,1])))

#Create a model profile, by convolving a Gaussian instrument profile with 
#a Gaussian line
gsig =  fwhm / 2.3548
gauss_prof = np.exp(-x_subsamp**2/2.0/( (line_width/km_sec_pix)**2 + gsig**2))

pixel_phase = np.linspace(-0.5,0.5,nphase)
fitpar = np.empty( (nphase,4) )
npix = 13
x_data = -7+np.arange(npix)
sig = np.ones(npix)/snr
for i in range(nphase):
    xc = pixel_phase[i] + x_data
    two_lines = np.interp(xc - line_sep/2.0/km_sec_pix,x_subsamp,line) + \
        line_rat*np.interp(xc + line_sep/2.0/km_sec_pix,x_subsamp,line)
    two_lines += np.random.normal(size=npix)/np.sqrt(km_sec_pix/2.0)/snr
    #Now fit two scaled gauss_prof functions to this.
    initp = np.array([-line_sep/2.0/km_sec_pix, 1.0, line_sep/2.0/km_sec_pix, 0.8])
    bestp = op.leastsq(two_line_resid,initp,args=(x_subsamp, gauss_prof, x_data, two_lines, sig))
    fitpar[i,:] = bestp[0]
    
plt.clf()
plt.plot(x_data, 1-two_lines*0.5,'o')
p1 = bestp[0][1]*np.interp(x_subsamp - bestp[0][0],x_subsamp,gauss_prof)*0.5
p2 = bestp[0][3]*np.interp(x_subsamp - bestp[0][2],x_subsamp,gauss_prof)*0.5
plt.plot(x_subsamp,1-p1)
plt.plot(x_subsamp,1-p2)
plt.plot(x_subsamp,1-(p2 + p1))
plt.axis([-4.5,4.5,0.6,1.05])
plt.xlabel('Pixels')
plt.ylabel('Flux')
    
print(km_sec_pix)
print("Sep: {0:5.3f} +/- {1:5.3f}".format(np.mean(fitpar[:,2]-fitpar[:,0])*km_sec_pix,np.std(fitpar[:,2]-fitpar[:,0])*km_sec_pix))
print("Rat: {0:5.3f} +/- {1:5.3f}".format(np.mean(fitpar[:,3]/fitpar[:,1]),np.std(fitpar[:,3]/fitpar[:,1])))