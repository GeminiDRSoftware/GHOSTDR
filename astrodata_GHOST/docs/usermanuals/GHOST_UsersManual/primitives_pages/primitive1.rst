.. primitive1:

.. rejectCosmicRays:

rejectCosmicRays
============================

Purpose
-------
This primitive will mask cosmic rays in a GHOST data frame by updating the
data mask.

.. warning:: The algorithm seems to work, but isn't managing to update the
             DQ plane. I'm working on the issue. -MCW

Inputs and Outputs
------------------
The following options can be passed to rejectCosmicRays via the RecipeSystem,
or overridden with a user configuration file, as marked:

+------------+---------+---------------------+-----------------+--------------------------------------+
| Parameter  |  Type   | Default             |     Override?   | Description                          |
+            +         +                     +---------+-------+                                      +
|            |         |                     | Recipe  |  User |                                      |
+============+=========+=====================+=========+=======+======================================+
| suffix     | str     | _cosmicRaysRejected | Y       | Y     | Suffix affixed to modified file.     |
+------------+---------+---------------------+---------+-------+--------------------------------------+
|subsampling | int     | 2                   | Y       | Y     | The number of pixels to subsample in |
|            |         |                     |         |       | each direction of each pixel.        |
+------------+---------+---------------------+---------+-------+--------------------------------------+
| sigma_lim  | float   | 5.0                 | Y       | Y     | Sigma clipping value.                |
+------------+---------+---------------------+---------+-------+--------------------------------------+
| f_lim      | float   | 6.0                 | Y       | Y     | Fine-structure clipping limit.       |
+------------+---------+---------------------+---------+-------+--------------------------------------+
| n_passes   | int     | 2                   | Y       | Y     | Number of iterations to perform.     |
+------------+---------+---------------------+---------+-------+--------------------------------------+

.. _ADS: https://ui.adsabs.harvard.edu/#abs/2001PASP..113.1420V/abstract

Algorithm
---------
This primitive is a first-principles implementation of the LACosmic algorithm.
NumPy array operations are used for maximum efficiency.
For a full-discussion of the algorithm, see van Dokkum 2001, PASP **113**,
1420-1427 (ADS_).
In summary (where the symbol :math:`\ast` denotes convolution):

1. Create a Laplacian, :math:`\mathcal{L}^+`, of the input data frame,
   :math:`I`.

   - Subsample the image by a factor of ``subsampling`` and apply the
     Laplacian transformation;
   - Set any negative sub-pixels to have value 0 instead;
   - Re-sample the image back to its original resolution.

2. Generate the 'sigma map', :math:`S`, of the frame.

   - Form a 'noise image' of the data,
     :math:`N=g^{-1}\sqrt{g(M_5\ast I)+\sigma_{rn}^2}`, where :math:`g` and
     :math:`\sigma_{rn}` are the data gain and read noise respectively.
     :math:`M_5` is a :math:`5\times 5` median filter used to smooth the data.
   - The sigma map is then generated as :math:`S=\mathcal{L}^+ /f_S N`, where
     :math:`f_S` is the ``subsampling`` factor used.
   - Smooth the sigma map to help remove sampling flux, resulting in
     :math:`S^\prime = S - (S\ast M_5)`.

3. Generate the 'fine-structure image' of the frame,
   :math:`\mathcal{F}=(M_3\ast I) - ((M_3\ast I) \ast M_7)`, where :math:`M_3`
   and :math:`M_7` are :math:`3\times 3` and :math:`7\times 7` median filters,

4. Mark pixels as cosmic rays (via the image quality plane) if:

   - :math:`S^\prime >` ``sigma_lim``, and
   - :math:`\mathcal{L}^+ /\mathcal{F} >` ``f_lim``.

5. For the purposes of the algorithm (i.e. **not** in the data to be returned),
   replace cosmic rays pixels with the median value of surrounding pixels.

6. Iterate over steps 1-5 until ``n_passes`` is reached, or the number of
   cosmic rays detected falls to 0.


Issues and Limitations
----------------------
Both strong spectral features and, in high-resolution mode, the ThXe calibration
output may be incorrectly flagged as cosmic rays. Therefore, this
primitive will trigger a call to
**need to add more here once this actually works**
