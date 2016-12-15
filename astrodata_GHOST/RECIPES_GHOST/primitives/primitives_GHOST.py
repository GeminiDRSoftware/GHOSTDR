from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors
from astrodata.utils.ConfigSpace import lookup_path
from gempy.gemini import gemini_tools as gt
from gempy.adlibrary.mosaicAD import MosaicAD
from gempy.gemini.gemMosaicFunction import gemini_mosaic_function
from astrodata_GHOST.ADCONFIG_GHOST.lookups import timestamp_keywords as ghost_stamps
from gempy.gemini.eti.gireduceparam import subtract_overscan_hardcoded_params
from pyraf import iraf
import numpy as np
import scipy
import functools

from astrodata_Gemini.RECIPES_Gemini.primitives.primitives_GMOS import GMOSPrimitives
from astrodata_Gemini.RECIPES_Gemini.primitives.primitives_stack import StackPrimitives

from astrodata_Gemini.ADCONFIG_Gemini.lookups import PolyfitDict

class GHOSTPrimitives(GMOSPrimitives):
    """
    Primary primitive set for GHOST.

    Attributes
    ----------
    astrotype : str
        Set to "GHOST"
    """

    astrotype = "GHOST"

    def init(self, rc):
        """
        GHOSTPrimitives init function.

        .. note:: This calls the GMOSPrimitives init function, changes
                  the timestamp_keys to GHOST ones, and returns the
                  ReductionContext.

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Returns
        -------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.
        """
        GMOSPrimitives.init(self, rc)
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)
        return rc

    def findApertures(self, rc):
        """
        Locate the apertures within a GHOST frame, and write out polyfit-
        compliant FITS tables to the calibrations system

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Returns
        -------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.
        """

        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "findAperture", "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["findAperture"]

        # Initialize the list of output AstroData objects
        # Note this will simply be the input list, with the timestamp_key
        # added
        adoutput_list = []

        # Loop over each file in the rc
        for ad in rc.get_inputs_as_astrodata():
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by findApertures"
                            % (ad.filename))
                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # This primitive should only be run on files of type GHOST_FLAT
            # which have been successfully prepared
            # Therefore, if these two types aren't in ad.types, skip the file
            if 'PREPARED' not in ad.types or 'GHOST_FLAT' not in ad.types:
                log.warning('findApertures is only run on prepared flats; '
                            'therefore, %s will not be used' % ad.filename)
                adoutput_list.append(ad)
                continue

            all_polyfit_dict = PolyfitDict.polyfit_dict
            # Work out the directory to get the Polyfit initial files from
            key = self._get_polyfit_key(ad)

            if key in all_polyfit_dict:
                # FIXME Hack currently required - can we not use lookup_path?
                # lookup_path assumes it's getting a single file, so will
                # prepend .py to any return that doesn't have a file type -
                # we need to strip this off to get a path
                poly_path = lookup_path(all_polyfit_dict[key]).replace('.py',
                                                                       '')
            else:
                # Don't know how to fit this file, so probably not necessary
                # Return this ad to the stream and move to the next
                log.warning('Not sure which initial model to use for %s; '
                            'skipping' % ad.filename)
                adoutput_list.append(ad)
                continue


            # Add the appropriate time stamps to the PHU
            gt.mark_history(
                adinput=ad, primname=self.myself(),
                keyword=timestamp_key)

            adoutput_list.append(ad)

        # Report the outputs to the RC
        rc.report_output(adoutput_list)

        yield rc

    def mosaicADdetectors(self, rc):
        """
        This primitive will mosaic the SCI frames of the input images, along
        with the VAR and DQ frames if they exist.

        :param tile: tile images instead of mosaic
        :type tile: Python boolean (True/False), default is False

        """
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "mosaicADdetectors", "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["mosaicADdetectors"]

        # Initialize the list of output AstroData objects
        adoutput_list = []

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():
            # Validate Data
            # if (ad.phu_get_key_value('GPREPARE')==None) and \
            #     (ad.phu_get_key_value('PREPARE')==None):
            #     raise Errors.InputError("%s must be prepared" % ad.filename)

            # Check whether the mosaicADdetectors primitive has been run
            # previously
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by mosaicADdetectors"
                            % (ad.filename))
                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # If the input AstroData object only has one extension, there is no
            # need to mosaic the detectors
            if ad.count_exts("SCI") == 1:
                log.stdinfo("No changes will be made to %s, since it " \
                            "contains only one extension" % (ad.filename))
                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # Get the necessary parameters from the RC
            tile = rc["tile"]

            log.stdinfo("Mosaicking %s ..." % ad.filename)
            log.stdinfo("MosaicAD: Using tile: %s ..." % tile)

            # trick mosaicing software into working for GHOST data
            ad.types.append('GSAOI')  # masquerade (briefly) as GSAOI
            ad.phu_set_key_value(
                'INSTRUME', '../../../astrodata_GHOST/ADCONFIG_GHOST/lookups')

            mo = MosaicAD(
                ad, mosaic_ad_function=gemini_mosaic_function,
                dq_planes=rc['dq_planes'])  # shift, rotate, and scale DQ plane
            adout = mo.as_astrodata(tile=tile)

            # undo hacks to trick mosaicing software
            ad.phu_set_key_value('INSTRUME', 'GHOST')
            adout.phu_set_key_value('INSTRUME', 'GHOST')
            ad.types.pop()
            ad.refresh_types()
            adout.refresh_types()

            # Add the appropriate time stamps to the PHU
            gt.mark_history(
                adinput=adout, primname=self.myself(), keyword=timestamp_key)

            # Change the filename
            adout.filename = gt.filename_updater(
                adinput=ad, suffix=rc["suffix"], strip=True)

            # Append the output AstroData object to the list
            # of output AstroData objects
            adoutput_list.append(adout)

        # Report the list of output AstroData objects to the reduction
        # context
        rc.report_output(adoutput_list)

        yield rc

    def rejectCosmicRays(self, rc):
        """
        Reject cosmic rays from GHOST data.

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
                                               function=np.median, footprint=fp,
                                               mode='constant',
                                               cval=np.nan)

            log.stdinfo('Doing CR removal for %s' % ad.filename)

            for i in range(ad['SCI'].count_exts()):
                amp = i + 1
                log.stdinfo('-----')
                log.stdinfo('AMP %d' % amp)
                log.stdinfo('-----')

                # Define an array that will hold the cosmic ray flagging
                # Note that we're deliberately not using the BPM at this stage,
                # otherwise the algorithm will start searching for cosmic rays
                # around pixels that have been flagged bad for another reason.
                cosmic_bpm = np.zeros_like(ad["SCI", amp].data,
                                           dtype=bool)

                # Start with a fresh copy of the data
                # Use numpy NaN to cover up any data detected bad so far
                # (i.e. BPM < 8)
                clean_data = np.copy(ad["SCI", amp].data)
                clean_data[ad['DQ', amp].data > 0] = np.nan

                no_passes = 0
                new_crs = 1
                new_cr_pix = None

                while new_crs > 0 and no_passes < rc['n_steps']:
                    no_passes += 1
                    curr_crs = np.count_nonzero(cosmic_bpm)
                    if curr_crs > 0 and new_cr_pix is not None:
                        # Median out the pixels already defined as cosmic rays
                        log.stdinfo('Pass %d: Median over previously '
                                    'found CR pix' % no_passes)

                        # One pass option - slow
                        # clean_data[new_cr_pix > 0] = median_replace(
                        #     clean_data)[new_cr_pix > 0]

                        # Loop option - faster for the number of CR (~ few k
                        # we expect for realistic data
                        inds = np.argwhere(new_cr_pix)
                        pad_data = np.pad(clean_data, 1, 'constant',
                                          constant_values=(np.nan, ))
                        # log.stdinfo('Padded array size: %s' % str(pad_data.shape))
                        # log.stdinfo(
                        #     'Data array size: %s' % str(clean_data.shape))
                        # log.stdinfo(
                        #     'CR array size: %s' % str(new_cr_pix.shape))
                        for ind in inds:
                            # log.stdinfo(str(ind))
                            # Using nanmedian stops nan values being considered
                            # in the ordering of median values
                            clean_data[zip(ind)] = np.nanmedian(
                                fp * pad_data[
                                     ind[0]:ind[0] + 3,
                                     ind[1]:ind[1] + 3
                                     ]
                            )

                    # Actually do the cosmic ray subtraction here
                    # ------
                    # STEP 1
                    # Construct a model for sky lines to subtract
                    # TODO: Add option for 'wave' keyword, which parametrizes
                    # an input wavelength solution function
                    # ------
                    log.stdinfo('Pass %d: Building sky model' % no_passes)
                    sky_model = scipy.ndimage.median_filter(clean_data,
                                                            size=[7, 1],
                                                            mode='constant',
                                                            cval=np.nan)
                    m5_model = scipy.ndimage.median_filter(clean_data,
                                                           size=[5, 5],
                                                           mode='constant',
                                                           cval=np.nan)
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
                    # log.stdinfo(
                    #     'data array size: %s' % str(data_shape))
                    subsampl_data = np.repeat(np.repeat(
                        ad["SCI", amp].data, subsampling, axis=1),
                        subsampling, axis=0
                    )
                    # log.stdinfo(
                    #     'subsampl_data array size: %s' % str(subsampl_data.shape))
                    # Convolve the subsampled data with the Laplacian kernel,
                    # trimming off the edges this introduces
                    # Bring any negative values up to 0
                    init_conv_data = scipy.signal.convolve2d(
                        subsampl_data, laplace_kernel)[1:-1, 1:-1]
                    init_conv_data[np.nonzero(init_conv_data <= 0.)] = 0.
                    # log.stdinfo(
                    #     'init_conv_data array size: %s' % str(init_conv_data.shape))
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
                    # log.stdinfo('conv_data array size: %s' % str(conv_data.shape))

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
                                                             size=[5, 5],
                                                             mode='constant',
                                                             cval=np.nan)
                    sig_detrend = sigmap - sig_smooth

                    # ------
                    # STEP 5
                    # Identify the potential cosmic rays
                    # ------
                    log.stdinfo('Pass %d: Flagging cosmic rays' % no_passes)
                    # Construct the fine-structure image
                    # (F, eqn 14 of van Dokkum)
                    m3 = scipy.ndimage.median_filter(subbed_data, size=[3, 3],
                                                     mode='constant',
                                                     cval=np.nan)
                    fine_struct = m3 - scipy.ndimage.median_filter(m3,
                                                                   size=[7, 7],
                                                                   mode=
                                                                   'constant',
                                                                   cval=np.nan)
                    # Pixels are flagged as being cosmic rays if:
                    # - The sig_detrend image (S') is > sigma_lim
                    # - The contrast between the Laplacian image (L+) and the
                    #   fine-structure image (F) is greater than f_lim
                    sigma_lim = rc['sigma_lim']
                    f_lim = rc['f_lim']
                    new_cr_pix = np.logical_and(sig_detrend > sigma_lim,
                                                (conv_data/fine_struct) >
                                                f_lim)
                    cosmic_bpm[new_cr_pix] = 1
                    new_crs = np.count_nonzero(cosmic_bpm) - curr_crs
                    log.stdinfo('Pass %d: Found %d CR pixels' %
                                (no_passes, new_crs, ))

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
        Stack GHOST frames using IRAF.

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


        Algorithms
        ----------
        This is a test of autodoc's abilities.

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


        Algorithms
        ----------
        This is a test of autodoc's abilities.

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

    def _get_polyfit_key(self, adinput=None):
        """
        Helper function - returns the path for finding initial polyfit models

        Parameters
        ----------
        adinput:
            Input AstroData object we wish to calibrate

        Returns
        -------
        path:
            The file path to the relevant polyfit models. There will be two
            files under each path, wavemod.fits and xmod.fits. The primitives
            calling _get_polyfit_key will need to open these files as
            required.
        """
        # The polyfit models are keyed by the instrument, and in the case of
        # GHOST, the arm, resolution mode, and date valid from.
        # Get the
        # instrument, the x binning and the y binning values using the
        # appropriate descriptors
        ad = adinput
        instrument = ad.instrument()
        detector_x_bin = ad.detector_x_bin()
        detector_y_bin = ad.detector_y_bin()
        if (instrument is None or detector_x_bin is None or
                    detector_y_bin is None):
            raise Errors.Error("Input parameters")

        key = '%s_%s_%s' % (instrument, detector_x_bin, detector_y_bin)

        if 'GHOST' in ad.types:
            # Need to get the arm and the mode as well to make the key
            arm = ad.arm()
            res_mode = ad.res_mode()
            key = '%s_%s_%s' % (key, arm, res_mode)

        return key




