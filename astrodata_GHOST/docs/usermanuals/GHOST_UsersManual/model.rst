.. model:

.. _GHOST_spectrograph_model:

***********************
Polynomial Model Method
***********************

Description of the Spectrograph Model principle
===============================================

The principle employed in the development of this pipeline relies heavily
on a polynomial principle for all the modelling of the spectrograph's
characteristics from the data. These include the location of the orders,
the wavelength scale and the reciprocal model that converts the sampled
slit from the slit viewer to its image on the spectrograph CCD.

The idea is that instead of a traditional empirical extraction where the
orders are "scanned" for their location and arc lines are detected blindly,
we form a model of the spectrograph and fit the parameters using the
flat fields and arcs to measure small changes in the spectrograph on a nightly basis.
Then, the extraction process becomes relatively trivial with knowledge
of where all the flux is and the wavelength scale is uniquely determined.
A simple advantage of this method is that the entire spectrograph is modelled as one,
instead of each individual order as a separate entity. Ultimately, measurements such as
radial velocity shifts can be determined using a single varying parameter,
as opposed to a combination of measured shifts in all orders.

The principle implemented is that of a sum of polynomials of polynomials.
This differs from an approach where each physical parameter of the spectrograph
is modelled individually and focusses on a series of coefficients that represent
various aspects of the CCD images. We thereby minimise the number of required
parameters that describe the data.


Mathematical principle
======================

The ``polyfit`` method uses files containing polynomial coefficients.

The functional form is:

        ..math::

            wave = p_0(m) + p_1(m)*y' + p_2(m)*y'^2 + ...)

        with :math: `y^{prime} = y - y_{\rm middle}`, and:

        .. math::

            p_0(m) = q_{00} + q_{01} * m' + q_{02} * m'^2 + ...

        with :math: `m^{prime} = m_{\rm ref}/m - 1`

        This means that the simplest spectrograph model should have:
        :math: `q_{00}` : central wavelength of order m_ref
        :math: `q_{01}` : central wavelength of order m_ref
        :math: `q_{10}` : central_wavelength/R_pix,
        with R_pix the resolving power / pixel.
        :math: `q_{11}` : central_wavelength/R_pix,
        with R_pix the resolving power / pixel.
        ... with everything else approximately zero.
