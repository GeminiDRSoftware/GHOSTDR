'''
/*
 *+
 * FUNCTION NAME: template<class T> int cosmic
 *
 * INVOCATION: cosmic(T* image, int nx, int ny, double exposed, float shieldcover,
 *                    float rate, bool mask, Image_Pixel P)
 *
 * PARAMETERS: (">" input, "!" modified, "<" output)
 * ! image - image buffer
 * > nx - number of pixels in x direction
 * > ny - number in y direction
 * > exposed - Number of seconds of integration for the image.
 * > shieldcover - Cutoff angle (in degrees) due to shield, assume shield is perfect
 * > rate - rate at which cosmic rays are incident (rays/s/cm/cm)
 * > mask - set if required to add area of effect mask
 * > P - size of each pixel in microns in 3 dimensions
 *
 * FUNCTION VALUE: Number of hits
 *
 * PURPOSE: generate simulated cosmic ray hits on a detector
 *
 * DESCRIPTION: Some assumptions:
 *  1.  All CR's liberate 100 (+/- 10) electrons per 0.1 micron travel.  10%
 *      of the cosmic rays are He nuclei, and thus have 4x the effect.
 *
 *  2.  The plane of the detector (+/- some extra) is blocked by a shield.
 *
 *  3.  The CR shield stops 100% of CRs, regardless of energy.
 *
 *  4.  A CR event in the detector has no lasting effect.
 *
 * EXTERNAL VARIABLES:
 *
 * PRIOR REQUIREMENTS:
 *
 * DEFICIENCIES:

 * ACKNOWLEDGMENT:
 *   Written by: Joel D. Offenberg, Raytheon ITSS
 *   Extra comments & code clean-up by Marc White, Research School of Astronomy
 *   & Astrophysics, The Australian National University
 *
 *-
 */
template<class T> int cosmic(T* image, int nx, int ny, double exposed, float shieldcover,
			      float rate, bool mask, Image_Pixel P)
{
  const float deg2rad = M_PI/180.;        // Number of radians per degree
  float n_rays;                           // # cosmic rays on detector
  int i;
  float xpos, ypos, zpos;                 // position of cosmic ray hit on detector
  float theta, cosphi,sinphi;             // direction angles (and cosine thereof) of cosmic ray velocity
  float dpath;                            // Number of microns per step.  In theory,
                                          // this is 1 by definition, except for the last
                                          // step, but instead we seem to be assuming
                                          // even not-quite-micron steps.
  float depath;                           // Number of electrons liberated per step.
  float dx, dy, dz;                       // Amount of motion per step
  float charge;                           // charge from a cosmic ray - use random number generator
  int NHe = 0;                            // count of He nucleii.
  float crmask[3][3] = {{1.06e-3,1.66e-2,1.06e-3}, // area of effect mask - the values are based on
			{1.66e-2,1.00,1.66e-2},    // Bernie Rauscher's cosmic ray measurements.
			{1.06e-3,1.66e-2,1.06e-3}};
  int xxpp,yypp;
  float crcont;
'''

import numpy as np
import math

# im_shape size of the image (2D)
# exposed  number of seconds of integration time
# shieldcover cutoff angle (in degrees) due to shield, measured as a zenith angle
# rate rate at which cosmic rays are incident (rays/s/cm/cm)
# use_mask True if the area effect mask is to be used
# pix_size size of pixels in microns (3D)


# returns image in electrons
def cosmic(im_shape, exposed, shieldcover, rate, use_mask, pix_size):
    # Input checking
    # Type-cast use_mask to make sure it's valid
    use_mask = bool(use_mask)

    image = np.zeros(im_shape)

    # Calculate number of rays based on rate for time exposed
    n_rays = im_shape[0] * pix_size[0] * 1e-4 * im_shape[1] * \
             pix_size[1] * 1e-4 * rate * exposed

    # Poisson number of cosmic rays
    n_rays = np.random.poisson(n_rays)
    # Make sure the number of cosmic rays is non-negative
    if n_rays < 0:
        n_rays = 0

    # These are the positions of the rays when they hit the detector
    pos = np.hstack((im_shape * np.random.random_sample((n_rays,2)),
                     np.zeros((n_rays,1))))
    # And their incoming angles, where theta is the azimuth angle and phi the zenith angle
    theta = 2.0 * math.pi * np.random.random_sample(n_rays)
    # We imagine that rays are more likely to come from angles closer to the zenith,
    # and we assume that the zenith is normal to the detector.
    sinphi = (np.random.random_sample(n_rays))**2
    cosphi = np.sqrt(1.0 - sinphi*sinphi)

    # Determine which rays are blocked by the shield
    hitargs = (sinphi > 0) * (sinphi < math.sin(math.radians(shieldcover)))
    pos = pos[hitargs]
    theta = theta[hitargs]
    sinphi = sinphi[hitargs]
    cosphi = cosphi[hitargs]

    # Number of electrons liberated per step
    depath = 100.0

    # Number of microns per step
    dpath = 0.1

    # Number of rays that hit the detector
    n_rays = len(pos)

    # 10% of the time the CR will be a He nucleus
    charge = np.where(np.random.random_sample(n_rays) > 0.9, 2.0, 1.0)

    # Area of effect mask - the values are based on
    # Bernie Rauscher's cosmic ray measurements.
    crmask = np.array([[1.06e-3, 1.66e-2, 1.06e-3],
                       [1.66e-2, 1.00,    1.66e-2],
                       [1.06e-3, 1.66e-2, 1.06e-3]])

    # The origin of the mask
    maskx0 = crmask.shape[0]/2
    masky0 = crmask.shape[1]/2

    # x, y, z elements per step
    dx = np.cos(theta) * dpath * sinphi
    dy = np.sin(theta) * dpath * sinphi
    dz = dpath * cosphi

    # Number of steps each CR takes through the detector
    # We ignore any fractional steps at the end of the path, trusting
    # that our step size is small enough to make the effect of these negligible
    nsteps = np.around(pix_size[2] / dz).astype(int)
    # the path each CR takes
    path = [
        (np.linspace(p[0], p[0]+n*x, n), np.linspace(p[1], p[1]+n*y, n))
        for p, x, y, n in zip(pos, dx, dy, nsteps)]

    # np.c_ (below) normally converts path into a 2D np.ndarray but sometimes,
    # in a low CR count regime, what results is a 3D np.ndarray (for reasons
    # still not understood); this 5 line hack forces the conversion to result
    # in the required 2D array
    shaper = np.ndarray(shape=(len(path), 2), dtype=np.ndarray)
    for i, p in enumerate(path):
        shaper[i][0] = p[0]
        shaper[i][1] = p[1]
    path = shaper

    path = np.c_[path, charge]

    # Add effect of CR at each step on its path
    for p in path:
        # Round the positions to integer pixels
        xs = np.around(p[0]).astype(int)
        ys = np.around(p[1]).astype(int)
        # Make sure they're in range
        args = (xs>=0) * (xs<im_shape[0]) * (ys>=0) * (ys<im_shape[1])
        xs = xs[args]
        ys = ys[args]
        if use_mask:
            for x,y in zip(xs,ys):
                crcont = np.random.poisson(depath*p[2]*p[2])
                # Add cosmic ray area-of-effect mask.
                for idx, mask in np.ndenumerate(crmask):
                    xx = x + idx[0] - maskx0
                    yy = y + idx[1] - masky0
                    if (xx >= 0) and (yy >= 0) and (xx < im_shape[0]) and (yy < im_shape[1]):
                        image[xx,yy] += (crcont * mask)[i]
        else:
            crcont = np.random.poisson(np.ones_like(xs)*depath*p[2]*p[2])
            if type(crcont) != np.ndarray:
                crcont = np.array([crcont])
            for i in range(len(crcont)):
                image[xs[i],ys[i]] += crcont[i]

    return image
