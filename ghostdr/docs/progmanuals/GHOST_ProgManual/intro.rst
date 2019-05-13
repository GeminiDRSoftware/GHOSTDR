.. intro:

******************************
Introduction to the Instrument
******************************

This section should be copied across from the "Overview of GHOST" section of
the User's Manual when that is finalized. The structure there is nearly
identical to the template here.

General Description
===================

.. note::
    This was copied directly from the Users' Manual on 10-05-19. Any updates
    made to the Users' Manual after that point need to be copied to here as
    well. -MCW

GHOST is a fibre-fed echelle spectrograph. It has a wide variety of observing
modes tuned to a range of science cases.

GHOST comprises two positioner arms, “IFU 1” and “IFU 2”. The layout of each
focal plane is shown below. Each comprises multiple micro-lens arrays
(black hexagons) feeding fiber bundles for targeting science fields at standard
(large hexagons in IFU 1 and 2) or high (small hexagons in IFU 1) spectral
resolution, peripheral guide fibers surrounding the science fields
(red hexagons), and dedicated sky fibers at a fixed offset position a few
arcseconds from the science fields (cyan hexagons).
Each positioner can access different halves of the 7.5 arcminute field
(as shown in Figure 4), with a 16” common region of overlap.

.. warning:: Find image

The instrument has two distinct observing modes. In 'standard' mode, up to two
targets can be observed simultaneously using the large fibre bundles on IFU1
and IFU2. In 'high resolution' mode, a single target can be observed using the
high-resolution bundle on IFU1. A ThXe calibration lamp can also be supplied
simultaneously with the observations. Sky fibre bundles are provided for
'standard' mode on IFU 2, and for 'high-resolution' mode on IFU 1.

+----------------------+-----------------------------------------+-------------------------+
| **Mode**             |          **Standard Resolution**        |    **High Resolution**  |
+======================+=========================================+=========================+
| **Spectral coverage**| 363-950nm, simultaneous                 | 363-950nm, simultaneous |
+----------------------+-----------------------------------------+-------------------------+
| **Spectral           | 50,000                                  | 75,000                  |
| resolution**         |                                         |                         |
+----------------------+-----------------------------------------+-------------------------+
| **Radial velocity    | 600 m/s                                 | 10 m/s                  |
| precision**          |                                         |                         |
+----------------------+-----------------------------------------+-------------------------+
| **Multiplexing**     | Dual targets/beam switching             | Single target           |
|                      | (82" minimum separation)                |                         |
+----------------------+-----------------------------------------+-------------------------+
| **Patrol field**     | 7.5' semicircle                         | 7.5'                    |
|                      | (16" overlap between IFUs)              |                         |
+----------------------+-----------------------------------------+-------------------------+
| **IFU plane size**   | 1.2"                                    | 1.2"                    |
+----------------------+-----------------------------------------+-------------------------+
| **IFU element        | 7 :math:`\times` 0.4"                   | 19 :math:`\times` 0.2"  |
| number and size**    |                                         |                         |
+----------------------+-----------------------------------------+-------------------------+
| **Sky fibres**       | 3 :math:`\times` 0.4"                   | 7 :math:`\times` 0.2"   |
|                      | (on IFU 1)                              | (on IFU 2)              |
+----------------------+-----------------------------------------+-------------------------+
| **Approx. limiting   | 17.4 - 17.8                             | 17.0 - 17.4             |
| Vega                 |                                         |                         |
| magnitude [1]_**     |                                         |                         |
+----------------------+-----------------------------------------+-------------------------+

.. [1] Can achieve S/N ratio of 30 in 1 hour at 450 nm.

The GHOST control software positions the IFUs as required during observations.
A mask in the slit injection unit blocks light from the unused fibre bundles.

Description of the Modes
========================

.. note::
    This was copied directly from the Users' Manual on 10-05-19. Any updates
    made to the Users' Manual after that point need to be copied to here as
    well. -MCW

he following table is an overview of the instrument observing modes available
to users, which are selecting via the Gemini Observing Tool (OT).

+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Configuration              | Focal Plane      | Resolution         | Detector          | Guiding [3]_      | Approx. B mag.     |
|                            |                  |                    | binning [2]_      |                   |                    |
+============================+==================+====================+===================+===================+====================+
| Two-target                 | IFU1 & IFU2      | 50k                | 1x2               | P2 & guide fibres | :math:`< 18`       |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Two-target (faint)         | IFU1 & IFU2      | 50k                | 1x4               | P2                | :math:`\geq 18`    |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Two-target (very faint)    | IFU1 & IFU2      | 43k                | 2x4               | P2                | :math:`\geq 18`    |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Beam-switching             | IFU1 & IFU2      | 50k                | 1x2               | P2 & guide fibres | :math:`< 18`       |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Beam-switching (faint)     | IFU1 & IFU2      | 50k                | 1x8               | P2                | :math:`\geq 18`    |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| Beam-switching (v. faint)  | IFU1 & IFU2      | 43k                | 2x8               | P2                | :math:`\geq 18`    |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| High-resolution            | IFU1             | 75k                | 1x1               | P2 & guide fibres | :math:`< 17`       |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| High-resolution (faint)    | IFU1             | 75k                | 1x8               | P2 & guide fibres | :math:`17-18`      |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+
| High-resolution (PRV)      | IFU1 & ThXe      | 75k                | 1x1               | P2 & guide fibres | :math:`<16`        |
+----------------------------+------------------+--------------------+-------------------+-------------------+--------------------+

.. [2] Reported as (binning in the spectral direction) x (binning in the slit
       direction)
.. [3] P2 corresponds to peripheral wave-front sensor (PWFS) guiding.

**Two-target modes**

In standard spectral resolution mode, two targets can be observed
simultaneously using the standard-resolution fibre bundles on IFU1 and IFU2,
and the associated sky bundle on IFU1. Targets are specified to the OT using
absolute astronomical coordinates; a guide star is also needed for the
peripheral wave-front sensor (PWFS) guiding. Care needs to be taken to
configure the instrument such that the PWFS does not vignette a science IFU.

For standard two-target modes, targets are expected to be bright enough that
the PSF edges can be used for guiding via the guide fibres attached to the
science bundles. For faint/very faint mode, guiding is by PWFS only. Guide
fibres can also be disabled in standard two-target mode if crowded fields
cause the guiding to be inaccurate.

In faint/very-faint mode, a larger detector binning is used to reduce the
impact of read noise.

**Beam-switching modes**

In regions of low target density (i.e. where there is a single target within
the GHOST field-of-view), the two standard-resolution IFUs may be beam-switched
to provide continuous target observation, whilst alternating each IFU between
the target and an offset sky position. This facilitates accurate sky
subtraction by differencing sequential frames, avoides the resampling of
bright sky lines or detector artefacts, and elimiates the effects of potential
flat fielding errors and differential fibre throughputs. This is particularly
useful for faint targets and the 'red' camera, where there are numerous
time-variable sky lines.

The Gemini OT will automatically set diametrically opposed offset conditions
for sky measurements, to allow beam-switching to be accomplished using
telescope motion alone. However, in the case that this is inappropriate (e.g.
crowded fields, or where the PWFS may vignette a science detector using
the default configuration), it is
possible to explicitly specify sky positions. This inflicts a time penalty, as
the IFU positions will need to be reconfigured.

In the faint and very-faint modes, larger detector binning is used, and
guiding via the GHOST fibres is disabled. Guide fibres may be used in the
standard beam-switching mode, although like the two-target mode, this can
be disabled if necessary.

**High-resolution modes**

High-resolution modes use the high-resolution science fibre bundle on IFU1. A
high-resolution sky fibre bundle is on IFU2, and can be positioned
independently of IFU1 for simultaneous sky observations. The use of a single
science field provides maximum flexibility for the positioning of IFUs so as
to avoid vignetting by the PWFS, and maximizes the patrol radius for selecting
PWFS guide stars. Spectral binning in the spectral direction is not used in this
mode, to fully sample the spectral PSF. A factor 2 binning along the slit is
optimal.

.. warning:: This factor 2 binning isn't reflected in the table!

The high-resolution science fibre bundle has six peripheral guide bundles, for
guiding using the extended PSF of bright targets. This can be disabled as
required, and is disabled by default in faint mode. Eight-pixel binning in the
slit direction is also used in faint mode.

For targets requiring the best possible wavelength calibration, a precision
radial velocity (PRV) mode is provided. A fibre agitator is used to reduce modal
noise introduced to the fibres by stress, strain or imperfections. A ThXe
calibration source is may also be fed into an additional high-resolution fibre
which is passed to the spectrograph for calibration simultaneous to
observations. This source is cycled on and off with a given duty cycle, giving
total counts within a given exposure time to be similar in magnitude to the
science fibres (and avoiding saturation).

**Spectropolarimetry mode**

.. note:: This mode is a desirable future upgrade.

In this mode, the two object probes are placed to one side of the field of
view under the spectropolarimetry module. A single star image is split into two
images in orthogonal polarization states (e.g., Stokes I+V and Stokes I-V),
with one probe detecting each polarization state. A standard acquisition
sequence is used to position each of the probes, and then multiple exposures
are taken with the polarization modulator in different states. For faint
sources, the two probes are beam switched fo that the sky fibres see the
difference in sky brightness at each output of the analyzer. In the
high-resolution mode, 50% of the light is lost, but observations are
otherwise identical.

Required Calibration and Associated Observations
================================================

Standard CCD Calibrators
------------------------

Biases and darks

Flat fields
-----------

Only fibre illuminated, not whole field - how is this then used?

Arc frames
----------

Need a before/after arc frame for each science observation, how these are
used etc.

Slit Calibrators
----------------

What kinds of slit calibrators are taken when, how used etc.

Important Instrument Characteristics and Issues
===============================================

This section cannot be completed until the instrument has been commissioned.

