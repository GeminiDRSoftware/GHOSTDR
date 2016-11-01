from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors
from gempy.gemini import gemini_tools as gt
from astrodata_GHOST.ADCONFIG_GHOST.lookups import timestamp_keywords as ghost_stamps
from gempy.gemini.eti.gireduceparam import subtract_overscan_hardcoded_params
from pyraf import iraf
import numpy as np
import scipy
import functools

from astrodata_Gemini.RECIPES_Gemini.primitives.primitives_GMOS import GMOSPrimitives
from astrodata_Gemini.RECIPES_Gemini.primitives.primitives_stack import StackPrimitives

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
        Reject cosmic rays from GHOST data.

        .. note:: This currently does not successfully flag anything in the DQ
                  plane. I'm looking into this. -MCW

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Yields
        ------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.

        """
        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "rejectCosmicRays",
                                 "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["rejectCosmicRays"]

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

            # Define the function for performing the median-replace of cosmic
            # ray pixels
            # Note that this is different from a straight median filter, as we
            # *don't* want to include the central pixel
            fp = [[1, 1, 1],
                  [1, 0, 1],
                  [1, 1, 1]]
            median_replace = functools.partial(scipy.ndimage.generic_filter,
                                               function=np.median, footprint=fp)

            log.stdinfo('Doing CR removal for %s' % ad.filename)

            for i in range(ad['SCI'].count_exts()):
                amp = i + 1

                # Define an array that will hold the cosmic ray flagging
                # Note that we're deliberately not using the BPM at this stage,
                # otherwise the algorithm will start searching for cosmic rays
                # around pixels that have been flagged bad for another reason.
                cosmic_bpm = np.zeros_like(ad["SCI", amp].data, dtype=bool)

                # Start with a fresh copy of the data
                clean_data = np.copy(ad["SCI", amp].data)

                no_passes = 0
                new_crs = 1
                while new_crs > 0 and no_passes < rc['n_passes']:
                    no_passes += 1
                    curr_crs = np.count_nonzero(cosmic_bpm)
                    # Median out the pixels already defined as cosmic rays
                    log.stdinfo('Pass %d: Wiping over previously '
                                'found bad pix' % no_passes)
                    if curr_crs > 0:
                        clean_data[cosmic_bpm] = median_replace(
                            clean_data)[cosmic_bpm]

                    # Actually do the cosmic ray subtraction here
                    # ------
                    # STEP 1
                    # Construct a model for sky lines to subtract
                    # TODO: Add option for 'wave' keyword, which parametrizes
                    # an input wavelength solution function
                    # ------
                    log.stdinfo('Pass %d: Building sky model' % no_passes)
                    sky_model = scipy.ndimage.median_filter(clean_data,
                                                            size=[7, 1])
                    m5_model = scipy.ndimage.median_filter(clean_data,
                                                           size=[5, 5])
                    subbed_data = clean_data - sky_model

                    # ------
                    # STEP 2
                    # Remove object spectra
                    # FIXME: Waiting on working find apertures routine
                    # ------

                    # ------
                    # STEP 3
                    # Compute 2nd-order Laplacian of input frame
                    # This is 'curly L' in van Dokkum 2001
                    # ------
                    # Subsample the data
                    log.stdinfo('Pass %d: Computing Laplacian' % no_passes)
                    subsampling = rc["subsampling"]
                    data_shape = ad["SCI", amp].data.shape
                    subsampl_data = np.repeat(np.repeat(
                        ad["SCI", amp].data, subsampling, axis=1),
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
                    log.stdinfo('Pass %d: Constructing sigma map' % no_passes)
                    gain = ad.gain()
                    read_noise = ad.read_noise()
                    noise = (1.0 / gain) * ((gain * m5_model +
                                             read_noise**2)**0.5)
                    noise_min = 0.00001
                    noise[np.nonzero(noise <= noise_min)] = noise_min
                    # div by 2 to correct convolution counting
                    sigmap = conv_data / (subsampling * noise)
                    # Remove large structure with a 5x5 median filter
                    # Equation (13) of van Dokkum 2001, generates S'
                    sig_smooth = scipy.ndimage.median_filter(sigmap,
                                                             size=[5, 5])
                    sig_detrend = sigmap - sig_smooth

                    # ------
                    # STEP 5
                    # Identify the potential cosmic rays
                    # ------
                    log.stdinfo('Pass %d: Flagging cosmic rays' % no_passes)
                    # Construct the fine-structure image
                    # (F, eqn 14 of van Dokkum)
                    m3 = scipy.ndimage.median_filter(subbed_data, size=[3, 3])
                    fine_struct = m3 - scipy.ndimage.median_filter(m3,
                                                                   size=[7, 7])
                    # Pixels are flagged as being cosmic rays if:
                    # - The sig_detrend image (S') is > sigma_lim
                    # - The contrast between the Laplacian image (L+) and the
                    #   fine-structure image (F) is greater than f_lim
                    sigma_lim = rc['sigma_lim']
                    f_lim = rc['f_lim']
                    cosmic_bpm[np.logical_and(sig_detrend > sigma_lim,
                                              (conv_data/fine_struct) > f_lim)] = 1
                    new_crs = np.count_nonzero(cosmic_bpm) - curr_crs
                    log.stdinfo('Found %d CR pixels in pass %d' % (new_crs,
                                                                   no_passes, ))

                # TODO: Determine whether to alter pix or alter BPM
                # For the moment, go with Mike Ireland's suggestion to require
                # a BPM update
                curr_bpm = ad["DQ", amp].data
                curr_bpm[cosmic_bpm] = 8
                ad['DQ', amp].data = curr_bpm

            # Append the ad to the output list
            adoutput_list.append(ad)

        # Finish
        rc.report_output(adoutput_list)

        yield rc

    def stackFrames(self, rc):
        """
        Reject cosmic rays from GHOST data.

        .. note:: This is a wrapped for the standard stackFrames primitive,
                  which increases the IRAF verbosity to 2.

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Yields
        ------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.

        """
        iraf.setVerbose(value=2)
        return StackPrimitives.stackFrames(self, rc)

    def standardizeHeaders(self, rc):
        """
        Standardize the headers of GHOST data to Gemini standards.

        .. note:: This function currently just called the equivalent
                  GMOS function (from GMOSPrimitives).

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Yields
        ------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.

        """
        #log = logutils.get_logger(__name__)
        #log.debug(gt.log_message("primitive", "standardizeHeaders", "starting"))
        #timestamp_key = self.timestamp_keys["standardizeHeaders"]
        rc.run('standardizeGeminiHeaders')
        yield rc

    def subtractOverscan(self, rc):
        """
        Subtract overscan from GHOST data.

        .. note:: This currently adds an extra subtract_overscan parameter
                  to send on to IRAF, then calls the equivalent GMOS
                  primitive (from GMOSPrimitives).

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Yields
        ------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.

        """
        iraf.setVerbose(value=2)
        # supply extra required param in gireduce's ETI wrapper; can also
        # override DETTYPE defaults this way (the alternative is to extend
        # gireduce to support our DETTYPE)
        subtract_overscan_hardcoded_params['order'] = 1
        return GMOSPrimitives.subtractOverscan(self, rc)

