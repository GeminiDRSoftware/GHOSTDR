.. extraction:

.. _GHOST_extraction_weights:

**********************
Extraction and Weights
**********************

.. warning:: This chapter needs significant expansion and further detail.

Weighted Extraction Background
==============================

In any spectrograph extraction routine, pixels in the direction orthogonal to the 
wavelength axis have to be averaged or summed in order to maximise the
signal-to-noise ratio.
In the case of multi-fiber spectroscopy, where there is some overlap between neighboring
fibers, the optimal algorithm is described in
Sharp and Birchall (2010, PASA, 27, 91;
`available via ADS <http://adsabs.harvard.edu/abs/2010PASA...27...91S>`_).

For GHOST, the situation is a little more complex, because the "object profile" is not
fixed as it is for the case of single-fiber spectroscopy. Each object has a profile
that comes from the slit viewing camera, which has to be scaled in the
spectrograph :math:`x`-direction (and also in flux).

Weighted Extraction Principle in 1 Dimension
============================================

Following the Sharp and Birchall notation (section 5.2 therein), we will
consider the profiles
for object and sky along a single order. These profiles, normalised to :math:`1`
along the :math:`x`-direction, are given by :math:`\phi_{ki}` for
object (or sky) :math:`k`, and :math:`x`-direction pixel :math:`i`.
From the variance plane (:math:`\sigma_i^2`) of the data, we obtain weights
for each pixel:

.. math::
    w_i = \frac{1}{\sigma_i^2},
    
with :math:`w_i=0` for bad pixels, including cosmic ray locations.

These weights are used to construct the following two matrices. First, the "cross-talk"
matrix between profiles, which should be diagonal except for the (very large!) 
cross-talk between sky and object:

.. math::
    C = c_{kj} = \Sigma_i \phi_{ki} \phi_{ji} w_i\textrm{,}

and the naive extraction matrix :

.. math::
    B = b_{kj} = \phi_{kj} w_j\textrm{.}
    
In the case of no cross-talk, this extraction matrix would left multiply the optimally 
weighted pixel values :math:`\mathbf{d}=D_i` to solve for the true flux :math:`\eta_k`:

.. math::
    \eta_k = \frac{\Sigma_{i} \phi_{ki} w_i D_i}{\Sigma_{i} \phi_{ki}^2 w_i} = \frac{B \cdot \mathbf{d}}{\Sigma_{i} \phi_{ki}^2 w_i} \textrm{.}

In the presence of cross-talk, the simple division on the right hand side of this equation
is not possible, and we aim to solve: 

.. math::
    C \cdot \mathbf{\eta} = B \cdot \mathbf{d},
    
which is solved by

.. math::
    \eta_k = \Sigma_i z_{ki} D_i
    
for an extraction weights matrix :math:`Z=z_{ki}` being computed by

.. math::
    Z = C^{-1} \cdot B\textrm{.}

This matrix also subtracts the sky as measured. This may not always be signal-to-noise
optimal. For example, if emission from the sky or object 2 is negligible, object 1 
is best extracted by ignoring those objects. For this reason, data are extracted in modes
corresponding to models with 1, 2 and 3 objects (including sky).

The extraction weights :math:`z_{kj}` 
are stored in an array that is :math:`n_x \times n_y` pixels
in size, where weights from every order are added together. This method of storing the
weights assumes no overlap in weights between orders - otherwise the weights from 
neighboring orders would be saved on the same pixels. This is a problem for 8x1 binning
in the shortest-wavelength orders, where there is overlap (at the level of up to 1 binned 
pixel) between neighboring orders. In this case, it is essential that these pixels are 
marked as `bad' in the reduction in both orders, to avoid cross-talk between orders. 

Flat Field Correction
=====================

There are two primary algorithms for flat field correction in spectroscopy: a long slit
flat field, and an extracted fiber profile flat field. GHOST is in an intermediate regime,
where the object profile varies (due to different fibers being illuminated), yet is
mostly stable. The primary flat field algorithm is similar to a long-slit flat field. 
During flat field processing, a model flat profile :math:`\phi_{k}` is constructed for 
every pixel :math:`k` for each of the red and blue cameras. This profile is normalised 
to 1.0 for every spectral direction pixel and every order. Note that we do not include 
a per-object index here, as the flat lamp can only illuminate the input fibers in the
same way as the sky. The median combined flat field is divided by the median for all
illuminated pixels, with no normalisation from order to order. This combined
flat field is defined as :math:`f_k`.

Prior to extraction, science data are corrected by the flat field as follows:

.. math::
    D_k = \frac{r_k \phi_{k}}{f_k}

where :math:`r_k` is the raw data for pixel :math:`k`.
Note that the Th/Xe simultaneous reference lamp is not flat-field corrected using this
algorithm, with those pixels treated separately. The Th/Xe lamp is only extracted in 1x1
binning mode for precision radial velocity observations.

Slit Tilt and PSF Variation Issues
==================================

One problem with the algorithm above is that it assumes that the two-dimensional 
point-spread-function (PSF) is the product of spectral and spatial point spread functions.
This is the case for single-fiber profiles imaged with spectrographs with typical
aberrations as used by Sharp and Birchall, producing a near-Gaussian profile. 
However, for GHOST, even individual
fiber profiles are not quite Gaussian. We avoid the most significant problems with this
while preserving the need to weight differently at high versus low signal-to-noise by
convolving the variance array both spatially and spectrally. In the current implementation,
the spectral weights are convolved by a Hanning window of ~7 pixels FWHM, and the spatial
variance is convolved by a (Hanning?) function of ~2 pixels FWHM.

Weighted Extraction Principle in 2D
===================================

For object profiles in two dimensions, we choose to extract as in one dimension, but
we interpolate linearly between the nearest two pixels in the wavelength direction only, 
so that the mean wavelength extracted for each pixel in the :math:`x`-direction corresponds
to the wavelength solution wavelength for the slit central pixel. This has the effect
of adding an additional dimension to the extraction weights variables:

.. math::
    \eta_k = \Sigma_l \Sigma_i z_{kil} a_{il} D_i\textrm{,}

where the sum over :math:`l` is over pixels in the dispersion direction, and the dispersion
interpolation weights are given by :math:`a_{il}`.
