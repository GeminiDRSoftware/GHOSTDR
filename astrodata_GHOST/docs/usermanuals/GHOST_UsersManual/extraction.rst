.. extraction:

.. _GHOST_extraction_weights:

**********************
Extraction and Weights
**********************

Weighted Extraction Background
==============================

In any spectrograph extraction routine, pixels in the direction orthogonal to the 
wavelength axis have to be averaged or summed in order to maximise the signal to noise.
In the case of multi-fiber spectroscopy where there is some overlap between neighboring
fibers, the optimal algorithm is described in Sharp and Birchall (2010, PASA, 27, 91).

For GHOST, the situation is a little more complex, because the "object profile" is not
fixed, as it is for the case of single fiber spectroscopy. Each object has a profile 
that comes from the slit viewing camera, which has to be scaled in the spectrograph `x`
direction and also in flux.

Weighted Extraction Principle in 1D
===================================

Following the Sharp and Birchall notation (section 5.2), we will consider the profiles 
for object and sky along a single order. These profiles, normalised to 1 along the x 
direction, are given by :math:`\phi_{ki}` for object (or sky) k and x-direction pixel i. 
From the variance plane of the data, we obtain weights for each pixel:

.. math::
    w_i = \frac{1}{\sigma_i^2},
    
with :math:`w_i=0` for bad pixels, including cosmic ray locations.

These weights are used to construct the following two matrices. First, the "cross-talk"
matrix between profiles, which should be diagonal except for the (very large!) 
cross-talk between sky and object:

.. math::
    C = c_{kj} = \Sigma_i \phi_{ki} \phi_{ji} w_i

and the naieve extraction matrix (which would be the extraction matrix in the case of 
no cross-talk):

.. math::
    B = b_{kj} = \phi_{kj} w_j

This is slightly different to the Sharp and Birchall notation, and for a data vector 
:math:`d=D_i`, their equation 11 becomes:

.. math::
    C \cdot \eta = B \cdot d,
    
which is solved by:

.. math::
    \eta_k = \Sigma_i z_{ki} D_i
    
for an extraction weights matrix :math:`Z=z_{ki}` being computed by:

.. math::
    Z = C^{-1} \cdot B

This matrix also subtracts the sky as measured. The extraction weights :math:`z_{kj}` 
are stored in an array that is nx :math:`\times` ny pixels 
in size, where weights from every order are added together. This only causes a potential
problem for 8x1 binning in the shortest-wavelength orders. For these orders, where there
is an overlap, it is important that those pixels are marked as `bad` in the reduction 
software.

Slit Tilt and PSF Variation Issues
==================================

One problem with the algorithm above is that it assumes that the two-dimensional 
point-spread-function is the product of spectral and spatial point spread functions.
This was the case for single fiber profiles imaged with spectrographs with typical
aberrations, producing a near-Gaussian profile. However, for GHOST, even individual
fiber profiles are not quite Gaussian. We avoid the most significant problems with this
while preserving the need to weight differently at high versus low signal-to-noise by
convolving the variance array.

Weighted Extraction Principle in 2D
===================================

For object profiles in 2-dimensions, we choose to extract as in 1-dimension, but where
we interpolate linearly between the nearest two pixels in the wavelength direction only, 
so that the mean wavelength extracted for each pixel in the `x` direction corresponds 
to the wavelength solution wavelength for the slit central pixel. This has the effect
of adding an additional dimension to the extraction weights variables:

.. math::
    \eta_k = \Sigma_l \Sigma_i z_{kil} a_{il} D_i

where the sum over l is over pixels in the dispersion direction, and the dispersion 
interpolation weights are given by :math:`a_{il}`