import math
import numpy as np
from skimage import transform, util

# pylint: disable=maybe-no-member, too-many-instance-attributes

'''
# Rotation of slit on camera
ROTANGLE = 90.0-9.04

# Point around which to rotate
ROT_CENTER = [900.0, 780.0]

SLITVIEW_PARAMETERS = {
    'std': {
        'central_pix': {
            'red': [781, 985],
            'blue': [771, 946]
            #'red': [77, 65],    #
            #'blue': [77, 156]   #
        },
        'extract_half_width': 3,
        'sky_pix_only_boundaries': {
            'red': [47, 63],
            'blue': [47, 63]
        },
        'object_boundaries': {
            'red': [[3, 46], [64, 107]],
            'blue': [[3, 46], [64, 107]]
        },
        'sky_pix_boundaries': {
            'red': [3, 107],
            'blue': [3, 107]
        },
    },
    'high': {
        'central_pix': {
            'red': [770, 857],
            'blue': [760, 819],
            #'red': [78, 95], #
            #'blue': [78, 4]  #
        },
        'extract_half_width': 2,
        'sky_pix_only_boundaries': {
            'red': [82, 106],
            'blue': [82, 106]
        },
        'object_boundaries': {
            'red': [[11, 81], [4, 9]],
            'blue': [[11, 81], [4, 9]]
        },
        'sky_pix_boundaries': {
            'red': [11, 106], 'blue': [11, 106]},
    }
}
'''

# A quick counterpoint to the existing floordiv (which is just a // b)
def ceildiv(a, b):
    return -(a // -b)

class SlitView(object):
    """
    A class containing tools common to processing the dark and bias corrected
    slit-viewer images.

    Parameters
    ----------
    slit_image: :obj:`numpy.ndarray`
        A single slit viewer image, which has been processed, cosmic-ray
        rejection etc.

    flat_image: :obj:`numpy.ndarray`
        A single slit viewer flat field image, which has been processed,
        cosmic-ray rejection etc.

    microns_pix: float (optional)
        Scale in microns in the slit plane for each pixel in the slit view-
        ing camera. The default value assumes 2x2 binning of the slit view-
        ing camera. Default is ``4.54*180/50*2``.

    mode: string
        ``'std'`` or ``'high'``: the observing mode. Default is ``'std'``.

    slit_length: float (optional)
        Physical slit length to be extracted in microns. Default is ``3600.``.
        
    reverse_profile: bool
        Do we reverse the profile? This is a sign convention issue between 
        the slit viewer and the CCD, to be determined through testing on 
        real data.
    """
    def __init__(self, slit_image, flat_image, slitvpars, microns_pix=4.54*180/50,
                 binning=2, mode='std', slit_length=3600.):
        self.binning = binning
        rota = slitvpars['rota']
        center = [slitvpars['rotyc'] // binning, slitvpars['rotxc'] // binning]
        self.central_pix = {
            'red': [slitvpars['center_y_red'] // binning,
                    slitvpars['center_x_red'] // binning],
            'blue': [slitvpars['center_y_blue'] // binning,
                     slitvpars['center_x_blue'] // binning]
        }
        self.sky_pix_only_boundaries = {
            'red': [slitvpars['skypix0'] // binning,
                    ceildiv(slitvpars['skypix1'], binning)],
            'blue': [slitvpars['skypix0'] // binning,
                     ceildiv(slitvpars['skypix1'], binning)]
        }
        self.object_boundaries = {
            'red': [[slitvpars['obj0pix0'] // binning,
                     ceildiv(slitvpars['obj0pix1'], binning)],
                    [slitvpars['obj1pix0'] // binning,
                     ceildiv(slitvpars['obj1pix1'], binning)]],
            'blue': [[slitvpars['obj0pix0'] // binning,
                     ceildiv(slitvpars['obj0pix1'], binning)],
                     [slitvpars['obj1pix0'] // binning,
                      ceildiv(slitvpars['obj1pix1'], binning)]],
        }
        if mode == 'std':
            self.sky_pix_boundaries = {
                'red': [slitvpars['obj0pix0'] // binning,
                        ceildiv(slitvpars['obj1pix1'], binning)],
                'blue': [slitvpars['obj0pix0'] // binning,
                        ceildiv(slitvpars['obj1pix1'], binning)]
            }
        else:
            self.sky_pix_boundaries = {
                'red': [slitvpars['obj0pix0'] // binning,
                        ceildiv(slitvpars['skypix1'], binning)],
                'blue': [slitvpars['obj0pix0'] // binning,
                        ceildiv(slitvpars['skypix1'], binning)]
            }
                            
        self.extract_half_width = ceildiv(slitvpars['ext_hw'], binning)
        if slit_image is None or rota == 0.0:
            self.slit_image = slit_image
        else:
            self.slit_image = transform.rotate(util.img_as_float64(slit_image), rota, center=center)
        if flat_image is None or rota == 0.0:
            self.flat_image = flat_image
        else:
            self.flat_image = transform.rotate(util.img_as_float64(flat_image), rota, center=center)
        self.mode = mode
        self.slit_length = slit_length
        self.microns_pix = microns_pix * binning
        self.reverse_profile = {'red': False, 'blue': True}
        
        # WARNING: These parameters below should be input from somewhere!!!
        # The central pixel in the y-direction (along-slit) defines the slit
        # profile offset, i.e. it interacts directly with the tramline fitting
        # and a change to one is a change to the other.
        # Co-ordinates are in standard python co-ordinates, i.e. y then x
        '''
        if mode in SLITVIEW_PARAMETERS.keys():
            for attr, value in SLITVIEW_PARAMETERS[mode].items():
                setattr(self, attr, value)
        else:
            raise ValueError("Invalid Mode: " + str(mode))
        '''

    def cutout(self, arm='red', use_flat=False):
        """
        Extract the 2-dimensional slit profile cutout.

        Parameters
        ----------
        arm: string, optional
            Either ``'red'`` or ``'blue'`` for GHOST. Default is ``'red'``.

        use_flat: bool, optional
            Cutout from the flat (True) or the slit frame (False). Default is
            False.

        Returns
        -------
        profile: :obj:`numpy.ndarray` (npix)
            The 2-dimensional slit profile cutout.
        """
        try:
            central_pix = self.central_pix[arm]
        except:
            raise ValueError("Invalid arm: '%s'" % arm)

        if use_flat:
            this_slit_image = self.flat_image
        else:
            this_slit_image = self.slit_image

        y_halfwidth = int(self.slit_length/self.microns_pix/2)
        return this_slit_image[
            central_pix[0]-y_halfwidth:central_pix[0]+y_halfwidth+1,
            central_pix[1]-self.extract_half_width:central_pix[1] +
            self.extract_half_width+1]

    def slit_profile(self, arm='red', return_centroid=False, use_flat=False,
                     denom_clamp=10, reverse_profile=None):
        """
        Extract the 1-dimensional slit profile.

        Parameters
        ----------
        arm: string, optional
            Either ``'red'`` or ``'blue'`` for GHOST. Default is ``'red'``.

        return_centroid: bool, optional
            Do we also return the pixel centroid of the slit? Default is False.
            
        use_flat: bool, optional
            Do we use the flat image? False for the object image.
            Default is False.
            
        denom_clamp: float, optional
            Denominator clamp - fluxes below this value are not used when
            computing the centroid. Defaults to ``10``.

        Returns
        -------
        profile: :obj:`numpy.ndarray` (npix)
            The summed 1-dimensional slit profile.
        """
        y_halfwidth = int(self.slit_length/self.microns_pix/2)
        cutout = self.cutout(arm, use_flat)
        if reverse_profile is None:
            reverse_profile = self.reverse_profile[arm]

        # Sum over the 2nd axis, i.e. the x-coordinate.
        profile = np.sum(cutout, axis=1)
        if reverse_profile:
            profile = profile[::-1]
        
        if return_centroid:
            xcoord = np.arange(
                -self.extract_half_width, self.extract_half_width+1)
            xcoord = np.tile(xcoord, 2*y_halfwidth+1).reshape(
                (2*y_halfwidth+1, 2*self.extract_half_width+1))
            centroid = np.sum(xcoord*cutout, 1)/np.maximum(profile, denom_clamp)
            return profile, centroid
        else:
            return profile

    def object_slit_profiles(self, arm='red', correct_for_sky=True,
                             used_objects=[0, 1],
                             append_sky=True, normalise_profiles=True):
        """
        Extract object slit profiles.

        Parameters
        ----------
        arm: string, optional
            Either ``'red'`` or ``'blue'`` for GHOST. Default is ``'red'``.

        correct_for_sky : bool, optional
            Should the slit profiles be corrected for sky? Defaults to True.

        append_sky : bool, optional
            Append the sky profile to the output ``profiles``? Defaults to True.

        normalise_profiles : bool, optional
            Should profiles be normalised? Defaults to True.
            
        used_objects: list of int, indices of used objects
            Denotes which objects should be extracted. Should be a list
            containing the ints 0, 1, or both, or None/the empty list
            to extract sky only.
            FIXME: Needs testing

        Returns
        -------
        profiles : list of :any:`numpy.ndarray`
            List of object slit profiles, as :any:`numpy.ndarray`.

        TODO: Figure out centroid array behaviour if needed.
        """
        # Input checking
        if used_objects is None:
            used_objects = []
        if len(used_objects) > 2:
            raise ValueError('used_objects must have length 1 or 2')
        used_objects = [int(_) for _ in used_objects]
        if not np.all([_ in [0, 1] for _ in used_objects]):
            raise ValueError('Only 0 and 1 may be in used_objects')
        if len(used_objects) != len(set(used_objects)):
            raise ValueError('Duplicate values are not allowed in '
                             'used_objects')

        # Find the slit profile.
        full_profile = self.slit_profile(arm=arm, reverse_profile=True)

        if correct_for_sky or append_sky:
            # Get the flat profile from the flat image.
            flat_profile = self.slit_profile(arm=arm, use_flat=True, reverse_profile=True)

        # WARNING: This is done in the extracted profile space. Is there any
        # benefit to doing this in pixel space? Maybe yes for the centroid.
        if correct_for_sky:
            flat_scaling = np.median(full_profile[
                self.sky_pix_only_boundaries[arm][0]:
                self.sky_pix_only_boundaries[arm][1] + 1
            ]) / np.median(flat_profile[
                self.sky_pix_only_boundaries[arm][0]:
                self.sky_pix_only_boundaries[arm][1] + 1
            ])
            full_profile -= flat_scaling*flat_profile

        # Extract the objects.
        profiles = []
        for boundary in [self.object_boundaries[arm][_] for _ in used_objects]:
            profiles.append(np.copy(full_profile))
            profiles[-1][:boundary[0]] = 0
            profiles[-1][boundary[1]+1:] = 0

        # Append the "sky" if needed (be aware that the top & bottom pixel
        # borders [the "edges" of the profile] contain some object flux)
        if append_sky:
            profiles.append(flat_profile)
            profiles[-1][:self.sky_pix_boundaries[arm][0]] = 0
            profiles[-1][self.sky_pix_boundaries[arm][1]+1:] = 0
        profiles = np.array(profiles)

        # Normalise profiles if requested (not needed if the total flux is what
        # you're after, e.g. for an mean exposure epoch calculation)
        if normalise_profiles:
            for prof in profiles:
                profsum = np.sum(prof)
                if math.isclose(profsum, 0.0, abs_tol=1e-9):
                    # FIXME what to do here?
                    raise ZeroDivisionError('Sum of profile is too close to zero')
                else:
                    prof /= np.sum(prof)

        if self.reverse_profile[arm]:
            return profiles
        else:
            return profiles[:,::-1]
