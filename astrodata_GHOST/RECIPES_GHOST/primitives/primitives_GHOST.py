from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors
from gempy.gemini import gemini_tools as gt
from astrodata_GHOST.ADCONFIG_GHOST.lookups import timestamp_keywords as ghost_stamps
from gempy.gemini.eti.gireduceparam import subtract_overscan_hardcoded_params
from pyraf import iraf
import numpy as np
import scipy

from primitives_GMOS import GMOSPrimitives
from primitives_stack import StackPrimitives

class GHOSTPrimitives(GMOSPrimitives):
    """
    Class containing all the GHOST primitives.  It inherits all the primitives
    from the GEMINIPrimitives class (which itself inherits a series of 
    primitives from RECIPE_Gemini/primitives.)
    """

    astrotype = "GHOST"

    def init(self, rc):
        GMOSPrimitives.init(self, rc)
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)
        return rc

    def rejectCosmicRays(self, rc):
        """
        Reject cosmic rays using a custom implementation of the LACosmic
        algorithm
        """
        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "rejectCosmicRays",
                                 "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["markCosmicRays"]

        # Initialise a list of output AstroData objects
        adoutput_list = []

        # Define the Laplacian and growth kernels for L.A.Cosmic
        laplace_kernel = np.array([
            [0.0, -1.0, 0.0],
            [-1.0, 4.0, -1.0],
            [0.0, -1.0, 0.0],
        ])
        growth_kernel = np.ones((3, 3), dtype=np.float64)

        # Loop over the input AstroData objects
        for ad in rc.get_inputs_as_astrodata():
            # Check if the rejectCosmicRays primitive has been run previously
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by rejectCosmicRays"
                            % ad.filename)

                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # Define an array that will hold the cosmic ray flagging
            # Note that we're deliberately not using the BPM at this stage,
            # otherwise the algorithm will start searching for cosmic rays
            # around pixels that have been flagged bad for another reason.
            cosmic_bpm = np.zeros_like(ad["SCI"].data, dtype=int)
            new_crs = 1

            while new_crs > 0:
                curr_crs = np.count_nonzero(cosmic_bpm)
                # Actually do the cosmic ray subtraction here
                # ------
                # STEP 1
                # Construct a model for sky lines to subtract
                # TODO: Add option for 'wave' keyword, which parametrizes
                # an input wavelength solution function
                # ------
                clean_data = np.copy(ad["SCI"].data)
                sky_model = scipy.ndimage.median_filter(clean_data, size=[7, 1])
                m5_model = scipy.ndimage.median_filter(clean_data, size=[5, 5])
                subbed_data = clean_data - sky_model

                # ------
                # STEP 2
                # Remove object spectra
                # TODO: Determine if this is necessary - PyWiFeS does not use it
                # ------

                # ------
                # STEP 3
                # Compute 2nd-order Laplacian of input frame
                # This is 'curly L' in van Dokkum 2001
                # ------
                # Subsample the data
                subsampling = rc["subsampling"]
                data_shape = ad["SCI"].data.shape()
                subsampl_data = np.repeat(np.repeat(
                    ad["SCI"].data, subsampling, axis=1),
                    subsampling, axis=0
                )
                # Convolve the subsampled data with the Laplacian kernel,
                # trimming off the edges this introduces
                # Bring any negative values up to 0
                init_conv_data = scipy.signal.convolve2d(
                    subsampl_data, laplace_kernel)[1:-1, 1:-1]
                init_conv_data[np.nonzero(init_conv_data <= 0.)] = 0.
                # Reverse the subsampling, returning the
                # correctly-convolved image
                conv_data = np.reshape(init_conv_data,
                                       (
                                           data_shape[0],
                                           init_conv_data.shape[0] //
                                           data_shape[0],
                                           data_shape[1],
                                           init_conv_data.shape[1] //
                                           data_shape[1],
                                       )).mean(axis=3).mean(axis=1)

                # ------
                # STEP 4
                # Construct noise model, and use it to generate the
                # 'sigma_map' S
                # This is the equivalent of equation (11) of van Dokkum 2001
                # ------
                gain = ad.gain()
                read_noise = ad.read_noise()
                noise = (1.0 / gain) * ((gain * m5_model + read_noise**2)**0.5)
                noise_min = 0.00001
                noise[numpy.nonzero(noise <= noise_min)] = noise_min
                # div by 2 to correct convolution counting
                # FIXME: Why is it divided by *two* times the noise?
                sigmap = conv_data / (2.0 * noise)
                # Remove large structure with a 5x5 median filter
                # Equation (13) of van Dokkum 2001, generates S'
                sig_smooth = scipy.ndimage.median_filter(sigmap, size=[5, 5])
                sig_detrend = sigmap - sig_smooth

                # ------
                # STEP 5
                # Identify the potential cosmic rays
                # ------
                # Construct the fine-structure image (F, eqn 14 of van Dokkum)
                m3 = scipy.ndimage.median_filter(subbed_data, size=[3, 3])
                fine_struct = m3 - scipy.ndimage.median_filter(m3, size=[7, 7])
                # Pixels are flagged as being cosmic rays if:
                # - The sig_detrend image (S') is > sigma_lim
                # - The contrast between the Laplacian image (L+) and the
                #   fine-structure image (F) is greater than f_lim
                sigma_lim = rc['sigma_lim']
                f_lim = rc['f_lim']
                cosmic_bpm[sig_detrend > sigma_lim and
                           (conv_data/fine_struct) > f_lim] = 1
                new_crs = np.count_nonzero(cosmic_bpm) - curr_crs

            # TODO: Determine whether to alter pix or alter BPM
            # For the moment, go with Mike Ireland's suggestion to require
            # a BPM update
            curr_bpm = ad["DQ"].data
            curr_bpm[cosmic_bpm] = True
            ad['DQ'].data = curr_bpm

            # Append the ad to the output list
            adoutput_list.append(ad)

        # Finish
        rc.report_output(adoutput_list)

        yield rc

    def stackFrames(self, rc):
        # Runs the standard stackFrames, but sets the iraf logging level higher
        iraf.setVerbose(value=2)
        return StackPrimitives.stackFrames(self, rc)

    def standardizeHeaders(self, rc):
        #log = logutils.get_logger(__name__)
        #log.debug(gt.log_message("primitive", "standardizeHeaders", "starting"))
        #timestamp_key = self.timestamp_keys["standardizeHeaders"]
        rc.run('standardizeGeminiHeaders')
        yield rc

    def subtractOverscan(self, rc):
        iraf.setVerbose(value=2)
        # supply extra required param in gireduce's ETI wrapper; can also
        # override DETTYPE defaults this way (the alternative is to extend
        # gireduce to support our DETTYPE)
        subtract_overscan_hardcoded_params['order'] = 1
        return GMOSPrimitives.subtractOverscan(self, rc)

