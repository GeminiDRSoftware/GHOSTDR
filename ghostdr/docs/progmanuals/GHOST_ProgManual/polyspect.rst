:tocdepth: 2

.. _polyspect:

***************
Polyspect Class
***************

.. note:: This is the first iteration of this documentation

This module contains the Polyspect class which contains the common tools for any
spectrograph. 

__init__
========

The class initialisation takes a series of inputs that define the spectrograph
orders and orientation

Most of the parameters are self-explanatory, but here is a list of those that
may not be:

+------------------------------+-----------------------------------------------+
| **Variable Name**            | **Purpose and meaning**                       |
+------------------------------+-----------------------------------------------+
| ``m_ref``, ``m_min`` and     | Reference, minimum and maximum order indeces  |
| ``m_max``                    | for the camera.                               |
+------------------------------+-----------------------------------------------+
| ``szx`` and ``szy``          | number of pixels in the x and y directions    |
+------------------------------+-----------------------------------------------+
| ``transpose``                | Whether the CCD orientation is transposed.    |
|                              | e.g. if the spectral direction is along y.    |
+------------------------------+-----------------------------------------------+
| ``x_map`` and                | The spatial and wavelength maps of the CCD.   |
| ``w_map``                    | These are populated by the spectral_format    |
+------------------------------+-----------------------------------------------+
| ``blaze`` and                | The blaze function and rotation matrices.     |
| ``matrices``                 | These are populated by the spectral_format too|
+------------------------------+-----------------------------------------------+


evaluate_poly
=============

This function is perhaps the key function that enables polynomial descriptions
of the various aspects of the spectrograph to work. In its simplest form, this
function converts the polynomial coefficients into values as a function of y
pixel and spectrograph order. The optional ``data`` input is only provided if
the desired structure of the output is not the default (orders, szy) evaluation
at all points.

This function is designed such that any set of polynomial coefficients can be
given and the evaluation will take place.

fit_resit
=========

This is the fit function for ``read_lines_and_fit`` and ``fit_to_x``.
A function of this type is needed by any least squares minimisation routine in
python, and this function takes the current model parameters ``params``, as well
as other needed values to evaluate the polynomial model, and returns the
weighted residuals to the fit.

The ``read_lines_and_fit`` function tries to minimise this residual and so does
the ``fit_to_x`` function. This was initially designed for the wavelength fit
only but later generalised for both wavelength and order position fitting.
The variable names are still named with the wavelength fitting mindset, but the
principle is the same. 

read_lines_and_fit
==================

This function reads an initial model parameter set and an array containing
spectral line information and fits to them.

This function takes an initial wavelength scale guess and a list of actual
arc lines and their wavelengths, sets up the least squares minimisation and fits
to the measured positions of the arc lines by the ``find_lines`` function in the
:doc:`extract` module.

spectral_format
===============

This function forms a spectrum from wavelength and a position polynomial model.

This is a simpler version of ``spectral_format_with_matrix`` (which replaces
this) that does not take the slit rotation into consideration.


adjust_x
========

This is an assisting function designed to slightly modify the initial position
fitting for any global shift in the spatial direction. The usefulness of this
function can be debated, but ultimately it won't create any problems.

From the function documentation:

Adjust the x pixel value based on an image and an initial
array from spectral_format().
This function performs a cross correlation between a pixel map
and a given 2D array image and calculates a shift along the spatial
direction. The result is used to inform an initial global shift for
the fitting procedure.
Only really designed for a single fiber flat, single fiber science
image or a convolution map.
This is a helper routine for :any:`fit_x_to_image`.


fit_x_to_image
==============

This is the main function that fits the tramline to the image. The initial
map has to be quite close and can be adjusted with the ``adjust_model`` function
and fitted manually and then fitted further.

This particular function sets up the data to be fitted, and then calls the
``fit_to_x`` to do the actual fit. Particularly, the ability and coding to
segment the image into smaller sections is done here, to limit the number of
parameters to fit and make the fitting faster. The default behaviour is to
subdivide the image by a factor of 8 in the spectral direction to reduce the
dimensions of the fit. This is fine because the final fit *must* be smooth
anyway.

This function also is responsible for taking the image to fit, typically the
result of the smoothed flat convolution, and find the location of the maxima
along the orders for fitting purposes. The process of cleaning up the flats is
done in the primitive itself, but here the weighting for each point is
determined based on the peak flux. This weighting is done to force the middle of
the orders to ensure the fainter edges are not contributing to the fit as much.

Optionally, the ``inspect`` parameter will display the initial order centre
maxima locations with point sizes following the weights, and then the result of
the fit.



fit_to_x
========

This is the helper function that takes the setup from ``fit_x_to_image`` and
actually fits the x model. This function also has a ``decrease_dim`` dimension
adjustment code, but is not used by default.

It's mostly self explanatory. 


spectral_format_with_matrix
===========================

This function is crucial to the workings of the pipeline, as it is designed to
be executed almost every time following the initialisation of the class as a way
to define the model variables for usage in the pipeline.

It takes the model parameters and defines the ``x_map``, ``w_map``, ``blaze``
and ``matrices`` arrays that describe every aspect of the spectrograph image for
extraction.

It makes use of the ``evaluate_poly`` function a few times to evaluate the
polynomials and then proceeds to use the models for slit rotation and spatial
and  spectral direction magnification to generate a (n_orders, y_pixel, 2, 2)
matrix that describes the slit rotation and magnification at each point along
an order.

The two orthogonal ``slit_microns_per_det_pix`` variables represent the physical
size of the full slit in detector pixels, which is scaled by the magnification
at all points. Each (2, 2) matrix contains 2 parameters that relate to the
magnification and are modified by rotation angle.

The mathematical principle behind this method is as follows. For every position
along an order we create a matrix :math:`A_{mat}` where we map the input angles
to output coordinates using the ``slit_microns_per_det_pix_x/y`` variables:

.. math::

   A_{mat} = \begin{bmatrix}
                1/slit\_microns\_per\_det\_pix\_y & 0 \\
                0 & 1/slit\_microns\_per\_det\_pix\_x
            \end{bmatrix}

We then compute an 'extra rotation matrix' :math:`R_{mat}` which maps the single
slit rotation value from the model for this position into a rotation matrix:

.. math::

   R_{mat} = \begin{bmatrix}
                cos(rot) & sin(rot) \\
                -sin(rot) & cos(rot)
            \end{bmatrix},

where :math:`rot` is the rotation of the slit in radians.

We then do a dot product of the two matrices and invert.

.. math::

   M = (R_{mat} . A_{mat})^{-1}

This is computationally complicated, but since the matrices are square and
:math:`A_{mat}` is diagonal, we can do this explicitly:

.. math::

   M = \begin{bmatrix}
   cos(rot) / slity & sin(rot) / slitx \\
   -sin(rot) / slity & cos(rot) / slitx
   \end{bmatrix} ^ {-1}

and

.. math::
   
   M = \begin{bmatrix}
   cos(rot) / slitx & -sin(rot) / slitx \\
   sin(rot) / slity & cos(rot) / slity
   \end{bmatrix}

Therefore, in the code we do this explicitly and obtain the result with simple
operations only. 

By default the function defines class variables and doesn't return anything.
Optionally, the actual arrays can be returned.


manual_model_adjust
===================

This function does quite complex things and is used mostly for visual inspection
and manual adjustment of the position and wavelength models.

It uses matplotlib slider widgets to adjust a polynomial model representation
overlaid on top of an image. A ``data`` array is provided
containing an image, and the model parameters needed to determine the desired
model. I.e., ``xparams`` are needed for both position and wavelength, and
``wparams`` are needed for the wavelength model.

A ``model`` variable is defined and needed to distinguish between a
representation of the order centre position or wavelength model.
A ``thar_spectrum`` array must be provided for the wavelength model.
The ``percentage_variation`` variable refers to the percentage of the parameter
values that should be allowed to range in the sliders. 

This function then goes on to show a window that contains the current
supplied model overlayed on top of the image provided, and a second window
containing a series of sliders which the user can use to vary the model
parameters in real time.

More documentation for this function may be added, but in principle the plotting
mechanism is commented and designed to be generic in terms of the number of
parameters. The actual workings of the plotting code is documented in the
matplotlib documentation.

The function returns the changed model parameters. 
