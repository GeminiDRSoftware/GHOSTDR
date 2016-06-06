************
Introduction
************

.. highlight:: python

PyGHOST is a pure python (i.e. no IRAF) analysis and simulation package for
GHOST data. The hope that is that it will become neat enough for easy 
incoporation into the Gemini recipe system.

 * reduction - the main module containing the reduction code.
 * simulation - the main module containing the simulation code.
 * doc - this documentation.
 * cal - calibration data

Basic Code Use
==============

Help for each of the modules inclues the simplest lines of code needed to run them.

Simulation Algorithms
=====================

The simulator is not a complete optical simulator - e.g. it doesn't directly deal with 
lenses and only propagates light according to first principles in collimated beams. The
two key equations in use are the grating equation, which in the plane orthogonal to the grating lines is:

.. math::

    n\lambda = d (\sin(\theta_i) - \sin(\theta_o)),

where :math: `n` is the order, :math: `\lambda` the wavelength, :math: `\theta_i` the input angle and :\math `theta_o` the output angle. In three dimensions, this becomes the following vector equation:

.. math::

    \mathbf{\hat{v}} \cdot \mathbf{\hat{s}} = \mathbf{\hat{u}} \cdot \mathbf{\hat{s}} + \frac{n\lambda}{d}

with :math:`\mathbf{s}` a unit vector perpendicular to the grating lines in the plane of the grating. Snell's law is also included in a similar way:

.. math::

    n_o \sin(\theta_o) &= n_i \sin(\theta_i) \\
    \mathbf{\hat{v}} &= \mathbf{\hat{n}} \cos(\theta_o) + \mathbf{\hat{p}} \sin(theta_i),
    
where :math: `\mathbf{\hat{n}}` is the surface normal, and `\mathbf{\hat{p}}` is a unit vector in both in the plane of the surface and in the plane defined by the input vector and surface normal.

Extraction Algorithms
=====================

The extraction algorithms are based on the mathematics in Sharp and Birchall (2009) 
(optimal extraction) and in Bolton and Schlegel (2009) (spectro-perfectionism). Neither
algorithm is used verbatim, because of the unique data analysis challenge of a long
fiber slit for GHOST that does not necessarily have adequate sampling. Instead, each
wavelength in each order (corresponding a single spectral-direction pixel in the center
of the slit image) has a nominal floating-point dispersion direction pixel coordinate
assigned to it. Then each row (i.e. dispersion direction) is multipled with a function
that is non-zero on up to 2 pixels, with a centroid equal to the pixel co-ordinate for that wavelength. Each row is then summed (i.e. the multiplication and sum is like a 
convolution) and the resulting column is treated exactly the same as in Sharp and Birchall
2009. 

Currently (July 2015), the sky and star extract very well, with equivalent widths and
resolution consistent between 1D and 2D extraction. However, the Th/Xe fiber gives the
same flux in 1D and 2D extraction for the flat lamp, but quite different fluxes whenever
there is a slit tilt for the arc lines. There is also some ringing in a Th/Xe fiber 
flat at extracted pixel
numbers higher than 3000 (where there is row tilt due to curvature), and a 
little deconvolution ringing evident in the Th/Xe 
2D extraction (also present, but to a lesser extent, in the 1D extraction). There 
may to be a ~0.1 pixel error in the extraction profile, but no more than that. So these
Th/Xe issues definitely seem to be an extraction artefact...
