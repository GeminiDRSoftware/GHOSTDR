import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# pylint: disable=maybe-no-member, too-many-instance-attributes


class SlitView(object):
    def __init__(self, slit_image, flat_image, microns_pix=4.54*180/50*2, mode='std',
                 slit_length=3600.):
        """
        A class containing tools common processing the dark and bias corrected
        slit-viewer images.

        Parameters
        ----------
        slit_image: numpy array
            A single slit viewer image, which has been processed, cosmic-ray
            rejection etc.

        microns_pix: float (optional)
            Scale in microns in the slit plane for each pixel in the slit view-
            ing camera. The default value assumes 2x2 binning of the slit view-
            ing camera.

        mode: string
            'std' or 'high': the observing mode.

        slit_length: float
            Physical slit length to be extracted in microns.

        !!!Lance - there need to be additional parameters here for you, stating
        the y-axis slit coordinate (or, if needed, y-axis pixels) corresponding
        to the boundaries of object 1 and 2 fluxes.
        """
        self.slit_image = slit_image
        self.mode = mode
        self.slit_length = slit_length
        self.microns_pix = microns_pix
        # WARNING: These parameters below should be input from somewhere!!!
        # The central pixel in the y-direction (along-slit) defines the slit
        # profile offset, i.e. it interacts directly with the tramline fitting
        # and a change to one is a change to the other.
        # Co-ordinates are in standard pythong co-ordinates, i.e. y then x
        if mode == 'std':
            self.central_pix = {'red': [84, 59], 'blue': [84, 150]}
            self.extract_half_width = 2
            #Boundaries for lower and upper pixels that contain *only* sky.
            #!!! WARNING: Change this.
            self.sky_pix_only_boundaries = {'red': [74,94], 'blue': [74,94]}
            #Boundaries for extracting the object. !!! WARNING: Change this.
            self.object_boundaries = {'red': [[30,75],[93,130]], 'blue': [[30,75],[93,130]]}
        elif mode == 'high':
            self.central_pix = {'red': [84, 95], 'blue': [84, 4]}
            self.extract_half_width = 3
            #Boundaries for lower and upper pixels that contain *only* sky.
            #!!! WARNING: Change this.
            self.sky_pix_only_boundaries = {'red': [34,54], 'blue': [34,54]}
            #!!! WARNING: Change this. The 2nd "object" from the point of view of
            #the extractor is the simultaneous Th/Xe. This could become an 
            #"object_type" parameter if we really cared.
            self.object_boundaries = {'red': [[53,75],[25,30]], 'blue': [[53,75],[25,30]]}
        else:
            raise UserWarning("Invalid Mode")

    def slit_profile(self, arm='red', return_centroid=False, use_flat=False):
        """Extract the 1-dimensional slit profile.

        Parameters
        ----------
        arm: string
            Either 'red' or 'blue' for GHOST.

        return_centroid: bool
            Do we also return the pixel centroid of the slit?

        Returns
        -------
        profile: numpy array (npix)
            The summed 1-dimensional slit profile.
        """
        try:
            central_pix = self.central_pix[arm]
        except:
            raise UserWarning("Invalid arm: " + arm)
        
        if use_flat:
            this_slit_image = self.flat_image
        else:
            this_slit_image = self.slit_image
            
        y_halfwidth = int(self.slit_length/self.microns_pix/2)
        cutout = this_slit_image[
            central_pix[0]-y_halfwidth:central_pix[0]+y_halfwidth+1,
            central_pix[1]-self.extract_half_width:central_pix[1]+
            self.extract_half_width+1]
        #Sum over the 2nd axis, i.e. the x-coordinate.
        profile = np.sum(cutout, axis=1)
        if return_centroid:
            xcoord = np.arange(
                -self.extract_half_width, self.extract_half_width+1)
            xcoord = np.tile(xcoord, 2*y_halfwidth+1).reshape(
                (2*y_halfwidth+1, 2*self.extract_half_width+1))
            centroid = np.sum(xcoord*cutout, 1)
            return profile, centroid
        else:
            return profile

    def object_slit_profiles(self, arm='red', correct_for_sky=True):
        """
        TODO: Figure out centroid array behaviour if needed.
        """
        object_boundaries = self.object_boundaries[arm]
        #Find the slit profile.
        full_profile = self.slit_profile(arm=arm)
        
        #WARNING: This is done in the extracted profile space. Is there any benefit to
        #doing this in pixel space? Maybe yes for the centroid.
        if correct_for_sky:
            #Get the flat profile from the flat image.
            flat_profile = self.slit_profile(arm=arm, use_flat=True)
            flat_scaling = np.median(full_profile[sky_pix_only_boundaries[0]:sky_pix_only_boundaries[1]+1])/
                np.median(flat_profile[sky_pix_only_boundaries[0]:sky_pix_only_boundaries[1]+1])
            full_profile -= flat_scaling*flat_profile
        
        #Extract the objects. 
        #WARNING: Dodgy code for now.
        profiles = []
        for boundary in object_boundaries:
            profiles.append(full_profile)
            profiles[-1][:boundary[0]]=0
            profiles[-1][boundary[1]+1:]=0
        
        return profiles