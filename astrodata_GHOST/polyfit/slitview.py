import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# pylint: disable=maybe-no-member, too-many-instance-attributes


class SlitView(object):
    def __init__(self, slit_image, microns_pix=4.54*180/50*2, mode='std',
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
        elif mode == 'high':
            self.central_pix = {'red': [84, 95], 'blue': [84, 4]}
            self.extract_half_width = 3
        else:
            raise UserWarning("Invalid Mode")

    def slit_profile(self, arm='red', return_centroid=False):
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
        y_halfwidth = int(self.slit_length/self.microns_pix/2)
        cutout = self.slit_image[
            central_pix[0]-y_halfwidth:central_pix[0]+y_halfwidth+1,
            central_pix[1]-self.extract_half_width:central_pix[1]+
            self.extract_half_width+1]
        profile = np.sum(cutout, 1)
        if return_centroid:
            xcoord = np.arange(
                -self.extract_half_width, self.extract_half_width+1)
            xcoord = np.tile(xcoord, 2*y_halfwidth+1).reshape(
                (2*y_halfwidth+1, 2*self.extract_half_width+1))
            centroid = np.sum(xcoord*cutout, 1)
            return profile, centroid
        else:
            return profile
