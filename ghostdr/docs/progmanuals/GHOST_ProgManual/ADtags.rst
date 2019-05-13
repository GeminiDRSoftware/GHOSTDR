.. ADtypes:

**************
AstroData Tags
**************

The following GHOST-specific tags are defined within :any:`ghost_instruments`.
Certain tags are only ever returned in pairs with another tags; such scenarios
are covered below.

``ghostdr`` AstroData Tags
==========================

Instrument identification
-------------------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``GHOST``          | Applies to all data taken | Automatically applied to    |
|                    | using the                 | any matched GHOST file      |
|                    | GHOST instrument          |                             |
+--------------------+---------------------------+-----------------------------+

.. note::
    All subsequent tags will only be applied to a file if the ``GHOST`` tag is
    applied.

Observation bundle
------------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``BUNDLE``         | A GHOST MEF bundle,       | Applied to any GHOST file   |
|                    | containing all data       | which does *not* match      |
|                    | relating to an            | any GHOST tag besides       |
|                    | observation               | ``GHOST``                   |
+--------------------+---------------------------+-----------------------------+

GCAL calibration observations
-----------------------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| (``CAL``,          | Bias calibration frame    | Header keyword              |
| ``BIAS``)          |                           | ``OBSTYPE == BIAS``         |
+--------------------+---------------------------+-----------------------------+
| (``CAL``,          | Dark calibration frame    | Header keyword              |
| ``DARK``)          |                           | ``OBSTYPE == DARK``         |
+--------------------+---------------------------+-----------------------------+
| (``CAL``,          | Flat calibration frame    | Header keyword              |
| ``FLAT``)          |                           | ``OBSTYPE == FLAT``         |
+--------------------+---------------------------+-----------------------------+
| (``CAL``,          | Arc calibration frame     | Header keyword              |
| ``ARC``)           |                           | ``OBSTYPE == ARC``          |
+--------------------+---------------------------+-----------------------------+

Ohter calibration observations
------------------------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``SKY``            | Sky calibration frame     | Header keyword              |
|                    |                           | ``OBSTYPE == SKY``          |
+--------------------+---------------------------+-----------------------------+

Resolution tags
---------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``HIGH``           | High-resolution mode      | Header keyword              |
|                    |                           | ``SMPNAME == HI_ONLY``      |
+--------------------+---------------------------+-----------------------------+
| ``STD``            | Standard resolution mode  | ``not`` of the above        |
+--------------------+---------------------------+-----------------------------+

Camera tags
-----------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| (``BLUE``,         | Blue spectrograph camera  | Header keyword              |
| ``SPECT``)         |                           | ``CAMERA == BLUE``          |
+--------------------+---------------------------+-----------------------------+
| (``RED``,          | Red spectrograph camera   | Header keyword              |
| ``SPECT``)         |                           | ``CAMERA == RED``           |
+--------------------+---------------------------+-----------------------------+
| (``SLITV``,        | Slit viewing camera       | Header keyword              |
| ``IMAGE``)         |                           | ``CCDNAME`` begins with     |
|                    |                           | ``'Sony-ICX674'``           |
+--------------------+---------------------------+-----------------------------+
Binning tags
------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``1x1`` or         | All data binned according | Header keyword ``CCDSUM``   |
| ``1x2`` or         | to the tag value          | is the same for all         |
| ``1x4`` or         |                           | extensions in the file      |
| ``1x8`` or         |                           |                             |
| ``2x4`` or         |                           |                             |
| ``2x8``            |                           |                             |
+--------------------+---------------------------+-----------------------------+
| ``NxN``            | Extensions within file    | ``not`` of the above        |
|                    | have different binnings   |                             |
+--------------------+---------------------------+-----------------------------+

Processing state tags
---------------------

+--------------------+---------------------------+-----------------------------+
| **AstroData Tag**  | **Observation Type**      | **Requirements**            |
+--------------------+---------------------------+-----------------------------+
| ``PROCESSED``      | Frame has been processed  | The presence of any of the  |
|                    |                           | following header keywords:  |
|                    |                           | ``PRSLITIM``, ``PRSLITBI``, |
|                    |                           | ``PRSLITDA``, ``PRSLITFL``, |
|                    |                           | ``PRWAVLFT``, ``PRPOLYFT``  |
+--------------------+---------------------------+-----------------------------+

Tag exclusion rules
-------------------

Certain tags will 'block' the application of other tags. These rules are as
follows:

+--------------------+---------------------------------------------------------+
| **AstroData Tag**  | **Subsequently blocked tags**                           |
+--------------------+---------------------------------------------------------+
| ``SLITV``          | ``SPECT``, ``BUNDLE``                                   |
+--------------------+---------------------------------------------------------+
| ``RED``, ``BLUE``  | ``BUNDLE``                                              |
+--------------------+---------------------------------------------------------+
