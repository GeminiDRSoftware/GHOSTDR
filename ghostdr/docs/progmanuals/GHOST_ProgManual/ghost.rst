:tocdepth: 2

.. _ghost:

**************
GhostArm Class
**************

.. note:: This is the first iteration of this documentation

This module contains the basic GhostArm class that defines the ghost
spectrograph arm settings. It initialises and inherits all attributes and
methods from Polyspect, which is the module that contains all spectrograph
generic functions.

__init__
=========

The class initialisation takes the arm, resolution mode and binning modes as
inputs and defines all needed attributes.

It starts by initialising the Polyspect class with the correct detector sizes
``szx`` and ``szy``, order numbers (``m_min``, ``m_max``) and whether the CCD is
transposed. Transposed in this case implies that the spectral direction is in
the x axis of the CCD image, which is the case for the GHOST data.

Most of the parameters are self-explanatory, but here is a list of those that
may not be:

+------------------------------+-----------------------------------------------+
| **Variable Name**            | **Purpose and meaning**                       |
+------------------------------+-----------------------------------------------+
| ``m_ref``, ``m_min`` and     | Reference, minimum and maximum order indeces  |
| ``m_max``                    | for the camera.                               |
+------------------------------+-----------------------------------------------+
| ``szx`` and ``szy``          | number of pixels in the x and y directions    |
|                              |                                               |
+------------------------------+-----------------------------------------------+
| ``nlenslets``                | Number of lenslets in the IFU                 |
+------------------------------+-----------------------------------------------+
| ``lenslet_high_size`` and    | Unused                                        |
| ``lenslet_std_size``         |                                               |
+------------------------------+-----------------------------------------------+


bin_data
========

This function is used to create a binned equivalent of a spectrograph image
array for the purposes of equivalent extraction. This function is now
implemented else where as the

.. class:: ghostdr.ghost.GHOST._rebin_ghost_ad

function and takes care of all the binning. It is essentially used because in
principle an unbinned flat and/or arc frame could be used to calibrate binned
data by artificially binning the calibration images.


slit_flat_convolve
==================

Function that takes a flat field image and a slit profile and convolves the two
in 2D. Returns result of convolution, which should be used for tramline fitting.

This function is key towards determining the centre of each order. Given the
potentially overlapping nature of fiber images in flat field frames, a
convolution method is employed with a sampled slit profile in which the center
of the order will, ideally, match the profile best and deveal the maximum of the
convolution.

A convolution map is then fed into a fitting function where the location of the
maxima in the map are found and a model is fit to determine a continuous
function describing the centre of the orders.

The function currently contains code related to convoling with a fixed synthetic
profile, which is not used. This is legacy code and sometimes used only for
testing purposes. The normal usage is to have the spatial scale model parameters
as inputs which determine the slit magnification as a function of order and
pixel along the orders. The flat convolution is done using only a smaller number
of orders (defaults to 3) and interpolated over the others but could in
principle be done with all orders considered. 

This function starts by taking a Fourier transform of the flat field
in the line::

  # Fourier transform the flat for convolution
  im_fft = np.fft.rfft(flat, axis=0)

It then proceeds to create a linear space for which orders to be evaluated::

  orders = np.linspace(self.m_min, self.m_max, num_conv).astype(int)

Then a convolution is done in 2D by interpolating the magnified slit profile
with the slit coordinates, normalising it and inverse Fourier transforming the
product between the flat transform and the shifted slit profile::

  #Create the slit model.
  mod_slit = np.interp(profilex*spat_scale[i], slit_coord, slit_profile)
  
  # Normalise the slit model and Fourier transform for convolution
  mod_slit /= np.sum(mod_slit)
  mod_slit_ft = np.fft.rfft(np.fft.fftshift(mod_slit))

  flat_conv_cube[j, :, i] = np.fft.irfft((im_fft[:, i] * mod_slit_ft)/num_conv)

Now, we have the convolution at ``num_conv`` locations, and the final result is
an interpolation between these.

