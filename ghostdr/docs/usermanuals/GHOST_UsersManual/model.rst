.. model:

.. _GHOST_spectrograph_model:

***********************
Polynomial Model Method
***********************

Description of the Spectrograph Model principle
===============================================

The GHOST data reduction software relies heavily
on a polynomial principle for all modelling of the spectrograph's
characteristics from the data. These include the location of the orders,
the wavelength scale, and the reciprocal model that converts the sampled
slit from the slit viewer to its image on the spectrograph CCD.

The idea is that, instead of a traditional empirical extraction (where the
orders are "scanned" for their location and arc lines are detected blindly),
we form a model of the spectrograph and fit the parameters using the
flat fields and arcs to measure small changes in the spectrograph on a nightly basis.
Then, the extraction process becomes relatively trivial with knowledge
of where all the flux is, and the wavelength scale is uniquely determined.
A major advantage of this method is that the entire spectrograph is modelled as one
cohesive unit,
instead of each individual order being treated as a separate entity.
Ultimately, this means measurements such as
radial velocity shifts can be determined using a single varying parameter,
as opposed to a combination of measured shifts in all orders.

The approach implemented is that of a sum of polynomials of polynomials.
This differs from an approach where each physical parameter of the spectrograph
is modelled individually, and focusses on a series of coefficients that represent
various aspects of the CCD images. We thereby minimise the number of required
parameters that describe the data.


Mathematical principle
======================

The :any:`polyfit` module produces and manipulates files containing polynomial
coefficients, where each row of the table are the coefficients for the
polynomials as a function of order. These are then combined as a function of
:math:`y`-position on the CCD chip. :math:`y` is defined as the CCD pixel
numbers in the spectral direction.

For mathematical convenience and correspondence with an easily-identified
reference, the polynomials
are evaluated with respect to:

- a reference order, :math:`m_{\rm ref}`, defaulting to whatever
  order number is in the middle of the range used for each arm; and,
- the middle pixel on the chip, :math:`y_{middle}`.

The functional form is:

.. math::

   F(p) = p_0(m) + p_1(m)*y' + p_2(m)*y'^2 + ...\textrm{,}

with :math:`y' = y - y_{middle}` , and:

        .. math::

	   p_0(m) = q_{00} + q_{01} * m' + q_{02} * m'^2 + ...\textrm{,}

with :math:`m' = m_{\rm ref}/m - 1`.

In this functional form, :math:`F(p)` is whatever aspect of the spectrograph
we wish to model. In the specific
example of GHOST, it will be the :math:`x`-position (defined in the spatial direction) in the first
instance, but this same method is then used for the wavelength scale, and all three aspects
of the slit image on the chip (spatial direction magnification scale, spectral direction
magnification scale and rotation). All of these aspects are expected to change
as a function of order, and position along the order.

This means that the simplest wavelength-scale spectrograph model should have:

 * :math:`q_{00}` : central wavelength of order :math:`m_{ref}`;
 * :math:`q_{01}` : central wavelength of order :math:`m_{ref}`;
 * :math:`q_{10}` : central_wavelength/:math:`R_{pix}`, with :math:`R_{pix}` the resolving power / pixel;
 * :math:`q_{11}` : central_wavelength/:math:`R_{pix}`, with :math:`R_{pix}` the resolving power / pixel;

with everything else approximately zero.

Note that the order of polynomials is left undefined. The code that handles these
parameters is identical and left generalised since each aspect (:math:`x` position, wavelength, etc)
may require a different number of variables to fully describe the spectrograph.

Description of model file contents
==================================

In the case of the :math:`x`-position, using default file ``xmod.fits``,
the contents of this file are as follows:

.. math::
  :label: e:matrix

  X_{mod} = \left[\begin{array}{ccccc}
     q_{24} & q_{23} & q_{22} & q_{21} & q_{20} \\
     q_{14} & q_{13} & q_{12} & q_{11} & q_{10} \\
     q_{04} & q_{03} & q_{02} & q_{01} & q_{00} \\
  \end{array}\right]

This non-standard way to define the variables within the files and imported array is related
to the way :any:`numpy`'s :any:`poly1d` function takes inputs,
with the highest-order coefficient first.

In the case of :math:`x`-position, the coefficients represent:

 * :math:`q_{00}` : :math:`x` position of the middle of the reference order;
 * :math:`q_{01}` : linear term coefficient for order spacing;
 * :math:`q_{02}` : quadratic term coefficient for order spacing;
 * :math:`q_{10}` : common rotation term for all orders;
 * :math:`q_{11}` : linear term coefficient for order rotation;
 * :math:`q_{12}` : quadratic term coefficient for order rotation;
 * :math:`q_{20}` : common curvature term for all orders;
 * :math:`q_{21}` : linear term coefficient for order curvature;
 * :math:`q_{22}` : quadratic term coefficient for order curvature;

with everything else approximately zero.
