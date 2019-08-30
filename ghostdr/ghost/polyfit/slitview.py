import numpy as np

# pylint: disable=maybe-no-member, too-many-instance-attributes


SLITVIEW_PARAMETERS = {
    'std': {
        'central_pix': {
            'red': [77, 65],
            'blue': [77, 156]
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
            'red': [78, 95],
            'blue': [78, 4]
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
    """
    def __init__(self, slit_image, flat_image, microns_pix=4.54*180/50*2,
                 mode='std', slit_length=3600.):
        self.slit_image = slit_image
        self.flat_image = flat_image
        self.mode = mode
        self.slit_length = slit_length
        self.microns_pix = microns_pix
        # WARNING: These parameters below should be input from somewhere!!!
        # The central pixel in the y-direction (along-slit) defines the slit
        # profile offset, i.e. it interacts directly with the tramline fitting
        # and a change to one is a change to the other.
        # Co-ordinates are in standard python co-ordinates, i.e. y then x
        if mode in SLITVIEW_PARAMETERS.keys():
            for attr, value in SLITVIEW_PARAMETERS[mode].items():
                setattr(self, attr, value)
        else:
            raise ValueError("Invalid Mode")

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
                     denom_clamp=10):
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

        # Sum over the 2nd axis, i.e. the x-coordinate.
        profile = np.sum(cutout, axis=1)
        if return_centroid:
            xcoord = np.arange(
                -self.extract_half_width, self.extract_half_width+1)
            xcoord = np.tile(xcoord, 2*y_halfwidth+1).reshape(
                (2*y_halfwidth+1, 2*self.extract_half_width+1))
            centroid = np.sum(xcoord*cutout, 1)/np.maximum(profile, denom_clamp)
            return profile, centroid
        else:
            return profile

    def object_slit_profiles(self, arm='red', correct_for_sky=True, used_objects=[0,1],
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
            
        used_objects: indices of used objects
            FIXME: Totally untested and handing off from Mike to Marc

        Returns
        -------
        profiles : list of :any:`numpy.ndarray`
            List of object slit profiles, as :any:`numpy.ndarray`.

        TODO: Figure out centroid array behaviour if needed.
        """
        # Find the slit profile.
        full_profile = self.slit_profile(arm=arm)

        if correct_for_sky or append_sky:
            # Get the flat profile from the flat image.
            flat_profile = self.slit_profile(arm=arm, use_flat=True)

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
        for boundary in [self.object_boundaries[arm][ix] for ix in used_objects]:
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
                prof /= np.sum(prof)

        return profiles
