#
#                                                                  gemini_python
#
#                                                      primitives_ghost_spect.py
# ------------------------------------------------------------------------------
import os
import numpy as np
from copy import deepcopy
import scipy
import scipy.signal as signal
import functools
from datetime import datetime, date, time, timedelta
import re
import astropy.coordinates as astrocoord
from astropy.time import Time
from astropy import units as u
from astropy import constants as const

import astrodata

from geminidr.gemini.lookups import DQ_definitions as DQ

from gempy.gemini import gemini_tools as gt
from gempy.mosaic.mosaicAD import MosaicAD

from .polyfit import GhostArm, Extractor, SlitView
from .polyfit.ghost import GhostArm

from .primitives_ghost import GHOST, filename_updater
from .parameters_ghost_spect import ParametersGHOSTSpect

from .lookups import polyfit_dict, line_list

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

GEMINI_SOUTH_LOC = astrocoord.EarthLocation.from_geodetic((-70, 44, 12.096),
                                                          (-30, 14, 26.700),
                                                          height=2722.,
                                                          ellipsoid='WGS84')

@parameter_override
class GHOSTSpect(GHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOST level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set(["GEMINI", "GHOST"])  # NOT SPECT because of bias/dark

    def __init__(self, adinputs, **kwargs):
        super(GHOSTSpect, self).__init__(adinputs, **kwargs)
        self.parameters = ParametersGHOSTSpect

    def addWavelengthSolution(self, adinputs=None, **params):
        """
        Apply the wavelength solution from an arc file (or set of arc files)
        to the data.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # No attempt to check if this primitive has already been run -
        # new arcs may be available which we wish to apply. Any old WAVL
        # extensions will simply be removed.

        # CJS: Heavily edited because of the new AD way
        # Get processed slits, slitFlats, and flats (for xmod)
        # slits and slitFlats may be provided as parameters
        arc_list = params.get("arc")
        # if arc_list is None:
        #     # CJS: This populates the calibrations cache (dictionary) with
        #     # "processed_slit" filenames for each input AD
        #     self.getProcessedArc(adinputs)
        #     # This then gets those filenames
        #     arc_list = [self._get_cal(ad, 'processed_arc')
        #                  for ad in adinputs]
        #     log.stdinfo(arc_list)

        # for ad, arcs in zip(
        #         *gt.make_lists(adinputs, arc_list, force_ad=True)):
        for i, ad in enumerate(adinputs):

            found_arcs = False

            if arc_list:
                try:
                    arc_before, arc_after = arc_list[i]
                    found_arcs = True
                except (TypeError, ValueError):
                    pass

            self.getProcessedArc(ad)
            if found_arcs == False:
                try:
                    arc_before, arc_after = self._get_cal(ad, 'processed_arc',)
                except (TypeError, ValueError):
                    # Triggers if only one arc, or more than two
                    arc_before = self._get_cal(ad, 'processed_arc',)[0]
                    arc_after = None

            log.stdinfo('Arcs for {}: {}, {}'.format(ad, arc_before, arc_after))

            # Stand up a GhostArm instance for this ad
            gs = GhostArm(arm=ad.arm(), mode=ad.res_mode(),
                          detector_x_bin=ad.detector_x_bin(),
                          detector_y_bin=ad.detector_y_bin())

            if arc_before is None:
                # arc = arc_after
                arc_after = astrodata.open(arc_after)
                wfit = gs.evaluate_poly(arc_after[0].WFIT)
            elif arc_after is None:
                # arc = arc_before
                arc_before = astrodata.open(arc_before)
                wfit = gs.evaluate_poly(arc_before[0].WFIT)
            else:
                # Need to weighted-average the wavelength fits from the arcs
                # Determine the weights (basically, the inverse time between
                # the observation and the arc)
                arc_after = astrodata.open(arc_after)
                arc_before = astrodata.open(arc_before)
                wfit_b = gs.evaluate_poly(arc_before[0].WFIT)
                wfit_a = gs.evaluate_poly(arc_after[0].WFIT)
                weight_b = np.abs((arc_before.ut_datetime() -
                                   ad.ut_datetime()).total_seconds())
                weight_a = np.abs((arc_after.ut_datetime() -
                                   ad.ut_datetime()).total_seconds())
                weight_a, weight_b = 1. / weight_a, 1 / weight_b
                log.stdinfo('Cominbing wavelength solutions with weights '
                            '%.3f, %.3f' %
                            (weight_a / (weight_a + weight_b),
                             weight_b / (weight_a + weight_b),
                             ))
                # Compute weighted mean fit
                wfit = wfit_a * weight_a + wfit_b * weight_b
                wfit /= (weight_a + weight_b)

            ad[0].WAVL = wfit

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def applyFlatBPM(self, adinputs=None, **params):
        """
        Find the flat relevant to the file(s) being processed, and merge the
        flat's BPM into the target file's.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        flat: str/None
            Name (full path) of the flatfield to use. If None, try:
        flatstream: str/None
            Name of the stream containing the flatfield as the first
            item in the stream. If None, the calibration service is used
        writeResult: bool
            Denotes whether or not to write out the result of profile
            extraction to disk. This is useful for both debugging, and data
            quality assurance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # No attempt to check if this primitive has already been run -
        # re-applying a flat BPM should have no adverse effects, and the
        # primitive simply skips if no flat is found.

        # CJS: extractProfile() contains comments explaining what's going on here
        flat_list = params["flat"]
        flat_stream = params["flat_stream"]
        if flat_list is None:
            if flat_stream is not None:
                flat_list = self.streams[flat_stream][0]
            else:
                self.getProcessedFlat(adinputs)
                flat_list = [self._get_cal(ad, 'processed_flat')
                            for ad in adinputs]

        for ad, flat in zip(*gt.make_lists(adinputs, flat_list, force_ad=True)):
            if flat is None:
                log.warning("No flat identified/provided for {} - "
                            "skipping".format(ad.filename))
                continue

            # Re-bin the flat if necessary
            # We only need the mask, but it's best to use the full rebin
            # helper function in case the mask rebin code needs to change
            if flat.detector_x_bin() != ad.detector_x_bin(
            ) or flat.detector_y_bin != ad.detector_y_bin():
                xb = ad.detector_x_bin()
                yb = ad.detector_y_bin()
                flat = self._rebin_ghost_ad(flat, xb, yb)
                # Re-name the flat so we don't blow away the old one on save
                flat_filename_orig = flat.filename
                flat.filename = filename_updater(flat,
                                                 suffix='_rebin%dx%d' %
                                                        (xb, yb,),
                                                 strip=True)
                flat.write(overwrite=True)

            # CJS: Edited here to require that the science and flat frames'
            # extensions are the same shape. The original code would no-op
            # with a warning for each pair that didn't, but I don't see how
            # this would happen in normal operations. The clip_auxiliary_data()
            # function in gemini_tools may be an option here.
            try:
                gt.check_inputs_match(adinput1=ad, adinput2=flat,
                                      check_filter=False)
            except ValueError:
                log.warning("Input mismatch between flat and {} - "
                            "skipping".format(ad.filename))
                continue

            for ext, flat_ext in zip(ad, flat):
                if ext.mask is None:
                    ext.mask = flat_ext.mask
                else:
                    ext.mask |= flat_ext.mask

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
            if params["write_result"]:
                ad.write(overwrite=True)

        return adinputs

    def barycentricCorrect(self, adinputs=None, **params):
        """
        Perform barycentric correction of the wavelength extension in the input
        files.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        correction_factor: float
            Barycentric correction factor to be applied. Defaults to None, at
            which point a computed value will be applied. The computed value
            is based on the recorded position of the Gemini South observatory.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by barycentricCorrect".
                            format(ad.filename))
                continue

            # TODO: Insert actual correction factor
            # Will the correction factor:
            # - Be allowed to be passed by the user?
            # - Be taken from the file header?
            # - Be computed here from the observing time of the image?
            if params.get('correction_factor') is None:
                cf = self._compute_barycentric_correction(ad, return_wavl=True)
            else:
                cf = params.get('correction_factor')

            # Multiply the wavelength scale by the correction factor
            for i, ext in enumerate(ad):
                log.stdinfo('Applying barycentric correction factor of '
                            '{} to ext {} of {}'.format(cf[i], i, ad.filename))
                ext.WAVL *= float(cf[i])

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def clipSigmaBPM(self, adinputs=None, **params):
        """
        Perform a sigma-clipping on the input data frame, such that any pixels
        outside the sigma threshold have their BPM value updated

        Parameters
        ----------
        sigma: float/None
            suffix to be added to output files
        bpm_value: int/None
            he integer value to be applied to the data BPM where the sigma
            threshold is exceeded. Defaults to 1 (which is the generic bad
            pixel flag). Note that the final output BPM is made using a
            bitwise_or operation.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        sigma = params["sigma"]
        bpm_value = params["bpm_value"]

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by clipSigmaBPM".
                            format(ad.filename))
                continue

            for ext in ad:
                extver = ext.hdr['EXTVER']
                if ext.mask is not None:
                    # TODO: Use astropy.stats.sigma_clip()
                    mean_data = np.mean(ext.data)
                    sigma_data = np.std(ext.data)
                    mask_map = (np.abs(ext.data-mean_data) > sigma*sigma_data)
                    if bpm_value:  # might call with None for diagnosis
                        ext.mask[mask_map] |= bpm_value

                    log.stdinfo('   {}:{}: nPixMasked: {:9d} / {:9d}'.format(
                        ad.filename, extver, np.sum(mask_map), ext.data.size))
                else:
                    log.warning('No DQ plane in {}:{}'.format(ad.filename,
                                                              extver))

            # Timestamp; DO NOT update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def darkCorrect(self, adinputs=None, **params):
        """
        GHOST-specific darkCorrect primitive

        This primitive, at it's core, simply copies the standard
        DRAGONS darkCorrect (part of :any:`Preprocess`). However, it has
        the ability to examine the binning mode of the requested dark,
        compare it to the adinput(s), and re-bin the dark to the
        correct format.

        To do this, this version of darkCorrect takes over the actual fetching
        of calibrations from subtractDark, manipulates the dark(s) as necessary,
        saves the updated dark to the present working directory, and then
        passes the updated list of dark frame(s) on to subtractDark.

        As a result, IOError will be raised if the adinputs do not
        all share the same binning mode.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        dark: str/list
            name(s) of the dark file(s) to be subtracted
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        sfx = params["suffix"]

        # Check if all the inputs have matching detector_x_bin and
        # detector_y_bin descriptors
        if not(all(
                [_.detector_x_bin() == adinputs[0].detector_x_bin() for
                 _ in adinputs])) or not(all(
            [_.detector_y_bin() == adinputs[0].detector_y_bin() for
             _ in adinputs]
        )):
            log.stdinfo('Detector x bins: %s' %
                        str([_.detector_x_bin() for _ in adinputs]))
            log.stdinfo('Detector y bins: %s' %
                        str([_.detector_y_bin() for _ in adinputs]))
            raise IOError('Your input list of files contains a mix of '
                          'different binning modes')

        # TODO: How to check for the primitive having already been applied
        # where the calibrations are bulk-fetched?
        # This is an issue, because if a dark has already been applied to the
        # file, it may no longer have an associated matching dark, which may
        # bomb the system.

        if params.get('dark', None):
            pass
        else:
            # All this line seems to do is check the valid darks can be found
            # for the adinputs
            self.getProcessedDark(adinputs)

        # Here we need to ape the part of subtractDark which creates the
        # dark_list, then re-bin as required, and send the updated dark_list
        # through to subtractDark
        # This is preferable to writing our own subtractDark, as it should
        # be stable against algorithm changes to dark subtraction
        dark_list = params["dark"] if params["dark"] else [
            self._get_cal(ad, 'processed_dark') for ad in adinputs]

        # We need to make sure we:
        # - Provide a dark AD object for each science frame;
        # - Do not unnecessarily re-bin the same dark to the same binning
        #   multiple times
        dark_list_out = []
        dark_processing_done = {}
        for ad, dark in zip(*gt.make_lists(adinputs, dark_list,
                                           force_ad=True)):
            if dark is None:
                if 'qa' in self.context:
                    log.warning("No changes will be made to {}, since no "
                                "dark was specified".format(ad.filename))
                    dark_list_out.append(None)
                    continue
                else:
                    raise IOError("No processed dark listed for {}".
                                  format(ad.filename))

            if dark.detector_x_bin() == ad.detector_x_bin() and \
                    dark.detector_y_bin() == ad.detector_y_bin():
                log.stdinfo('Binning for %s already matches input file' %
                            dark.filename)
                dark_list_out.append(dark.filename)
            else:
                xb = ad.detector_x_bin()
                yb = ad.detector_y_bin()
                dark = self._rebin_ghost_ad(dark, xb, yb)
                # Re-name the dark so we don't blow away the old one on save
                dark_filename_orig = dark.filename
                dark.filename = filename_updater(dark,
                                                 suffix='_rebin%dx%d' %
                                                        (xb, yb, ),
                                                 strip=True)
                dark.write(overwrite=True)
                dark_processing_done[
                    (dark_filename_orig, xb, yb)] = dark.filename
                dark_list_out.append(dark.filename)
                log.stdinfo('Wrote out re-binned dark %s' % dark.filename)

            # Check the inputs have matching binning, and shapes
            # Copied from standard darkCorrect (primitives_preprocess)
            # TODO: Check exposure time?
            try:
                gt.check_inputs_match(ad, dark, check_filter=False)
            except ValueError:
                # Else try to extract a matching region from the dark
                log.warning('AD inputs did not match - attempting to clip dark')
                dark = gt.clip_auxiliary_data(ad, aux=dark, aux_type="cal")

                # Check again, but allow it to fail if they still don't match
                gt.check_inputs_match(ad, dark, check_filter=False)

            log.stdinfo("Subtracting the dark ({}) from the input "
                        "AstroData object {}".
                        format(dark.filename, ad.filename))
            ad.subtract(dark)

            # Record dark used, timestamp, and update filename
            ad.phu.set('DARKIM', dark.filename, self.keyword_comments["DARKIM"])
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs

    def extractProfile(self, adinputs=None, **params):
        """
        Extract the object profile from a slit or flat image.
        
        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        slit: str/None
            Name of the (processed & stacked) slit image to use for extraction
            of the profile. If not provided/set to None, the primitive will
            attempt to pull a processed slit image from the calibrations
            database (or, if specified, the --user_cal processed_slit
            command-line option)
        slitflat: str/None
            Name of the (processed) slit flat image to use for extraction
            of the profile. If not provided, set to None, the RecipeSystem
            will attempt to pull a slit flat from the calibrations system (or,
            if specified, the --user_cal processed_slitflat command-line
            option)
        writeResult: bool
            Denotes whether or not to write out the result of profile
            extraction to disk. This is useful for both debugging, and data
            quality assurance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Make no attempt to check if primitive has already been run - may
        # have new calibrators we wish to apply.

        # CJS: Heavily edited because of the new AD way
        # Get processed slits, slitFlats, and flats (for xmod)
        # slits and slitFlats may be provided as parameters
        slit_list = params["slit"]
        if slit_list is None:
            # CJS: This populates the calibrations cache (dictionary) with
            # "processed_slit" filenames for each input AD
            self.getProcessedSlit(adinputs)
            # This then gets those filenames
            slit_list = [self._get_cal(ad, 'processed_slit')
                         for ad in adinputs]

        slitflat_list = params["slitflat"]
        if slitflat_list is None:
            self.getProcessedSlitFlat(adinputs)
            slitflat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        self.getProcessedFlat(adinputs)
        flat_list = [self._get_cal(ad, 'processed_flat') for ad in adinputs]

        # TODO: Have gt.make_lists handle multiple auxiliary lists?
        # CJS: Here we call gt.make_lists. This has only been designed to work
        # with one auxiliary list at present, hence the three calls. This
        # produces two lists of AD objects the same length, one of the input
        # ADs and one of the auxiliary files, from the list
        # of filenames (or single passed parameter). Importantly, if multiple
        # auxiliary frames are the same, then the file is opened only once and
        # the reference to this AD is re-used, saving speed and memory.
        _, slit_list = gt.make_lists(adinputs, slit_list, force_ad=True)
        _, slitflat_list = gt.make_lists(adinputs, slitflat_list, force_ad=True)
        _, flat_list = gt.make_lists(adinputs, flat_list, force_ad=True)

        for ad, slit, slitflat, flat in zip(adinputs, slit_list,
                                        slitflat_list, flat_list):
            log.info(ad.info())

            # CJS: failure to find a suitable auxiliary file (either because
            # there's no calibration, or it's missing) places a None in the
            # list, allowing a graceful continuation.
            if slit is None or slitflat is None or flat is None:
                log.warning("Unable to find calibrations for {}; "
                            "skipping".format(ad.filename))
                continue

            # CJS: Changed to log.debug() and changed the output
            log.debug("Slit parameters: ")
            log.debug("   processed_list: {}".format(slit.filename))
            log.debug("   processed_flat: {}".format(slitflat.filename))

            res_mode = ad.res_mode()
            arm = GhostArm(arm = ad.arm(), mode = res_mode,
                           detector_x_bin = ad.detector_x_bin(),
                           detector_y_bin = ad.detector_y_bin())

            # CJS: Heavy refactor. Return the filename for each calibration
            # type. Eliminates requirement that everything be updated
            # simultaneously.
            #key = self._get_polyfit_key(ad)
            #log.stdinfo("Polyfit key selected: {}".format(key))
            try:
                poly_wave = self._get_polyfit_filename(ad, 'wavemod')
                poly_spat = self._get_polyfit_filename(ad, 'spatmod')
                poly_spec = self._get_polyfit_filename(ad, 'specmod')
                poly_rot = self._get_polyfit_filename(ad, 'rotmod')
                wpars = astrodata.open(poly_wave)
                spatpars = astrodata.open(poly_spat)
                specpars = astrodata.open(poly_spec)
                rotpars = astrodata.open(poly_rot)
            except IOError:
                log.warning("Cannot open required initial model files for {};"
                            " skipping".format(ad.filename))
                continue

            arm.spectral_format_with_matrix(flat[0].XMOD, wpars[0].data,
                        spatpars[0].data, specpars[0].data, rotpars[0].data)
            sview = SlitView(slit[0].data, slitflat[0].data, mode=res_mode)

            extractor = Extractor(arm, sview, badpixmask=ad[0].mask,
                                  vararray=ad[0].variance)
            # CJS: Makes it clearer that you're throwing the first two
            # returned objects away (get replaced in the two_d_extract call)
            _, _, extracted_weights = extractor.one_d_extract(ad[0].data)
            extracted_flux, extracted_var = extractor.two_d_extract(ad[0].data,
                                    extraction_weights=extracted_weights)

            # CJS: Since you don't use the input AD any more, I'm going to
            # modify it in place, in line with your comment that you're
            # considering this.
            ad[0].reset(extracted_flux, mask=None, variance=extracted_var)
            ad[0].WGT = extracted_weights

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
            if params["write_result"]:
                ad.write(overwrite=True)

        return adinputs

    def interpolateAndCombine(self, adinputs=None, **params):
        """
        Perform a sigma-clipping on the input data frame, such that any pixels
        outside the sigma threshold have their BPM value updated
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by interpolateAndCombine".
                            format(ad.filename))
                continue

            # Timestamp; DO NOT update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def findApertures(self, adinputs=None, **params):
        """
        Locate the apertures within a GHOST frame, and insert a :any:`polyfit`
        model into a new extension on each data frame
        
        Parameters
        ----------
        slitflat: str/AD/None
            slit flat to use; if None, the calibration system is invoked
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Make no attempt to check if primitive has already been run - may
        # have new calibrators we wish to apply.

        # CJS: See comment in extractProfile() for handling of calibrations
        flat_list = params["slitflat"]
        if flat_list is None:
            self.getProcessedSlitFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        for ad, slit_flat in zip(*gt.make_lists(adinputs, flat_list, force_ad=True)):
            if not {'PREPARED', 'GHOST', 'FLAT'}.issubset(ad.tags):
                log.warning("findApertures is only run on prepared flats: "
                            "{} will not be processed".format(ad.filename))
                continue

            try:
                poly_xmod = self._get_polyfit_filename(ad, 'xmod')
                poly_spat = self._get_polyfit_filename(ad, 'spatmod')
                xpars = astrodata.open(poly_xmod)
                spatpars = astrodata.open(poly_spat)
            except IOError:
                log.warning("Cannot open required initial model files for {};"
                            " skipping".format(ad.filename))
                continue

            arm = ad.arm()
            res_mode = ad.res_mode()
            ghost_arm = GhostArm(arm=arm, mode=res_mode)

            # Create an initial model of the spectrograph
            xx, wave, blaze = ghost_arm.spectral_format(xparams=xpars[0].data)

            slitview = SlitView(slit_flat[0].data, slit_flat[0].data,
                                mode=res_mode)

            # This is an attempt to remove the worse cosmic rays
            # in the hope that the convolution is not affected by them.
            # Start by performing a median filter
            medfilt = signal.medfilt2d(ad[0].data, (5,5))
            # Now find which pixels have a percentage difference larger than
            # a defined value between the data and median filter, and replace
            # those in the data with the median filter values. Also, only
            # replace values above the data average, so as not to replace low
            # S/N values at the edges.
            data = ad[0].data.copy()
            condit = np.where(np.abs((medfilt - data)/(medfilt+1)) > 200) and\
                     np.where(data > np.average(data))
            data[condit] = medfilt[condit]

            # Convolve the flat field with the slit profile
            flat_conv = ghost_arm.slit_flat_convolve(data,
                slit_profile=slitview.slit_profile(arm=arm),
                spatpars=spatpars[0].data, microns_pix=slitview.microns_pix,
                xpars=xpars[0].data)

            flat_conv = signal.medfilt2d(flat_conv, (5,5))

            # Fit the initial model to the data being considered
            fitted_params = ghost_arm.fit_x_to_image(flat_conv,
                                                     xparams=xpars[0].data,
                                                     decrease_dim=8,
                                                     inspect=False)

            # CJS: Append the XMOD as an extension. It will inherit the
            # header from the science plane (including irrelevant/wrong
            # keywords like DATASEC) but that's not really a big deal.
            # (The header can be modified/started afresh if needed.)
            ad[0].XMOD = fitted_params

            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def fitWavelength(self, adinputs=None, **params):
        """
        Fit wavelength solution to a GHOST ARC frame
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Make no attempt to check if primitive has already been run - may
        # have new calibrators we wish to apply.

        self.getProcessedFlat(adinputs)
        flat_list = [self._get_cal(ad, 'processed_flat') for ad in adinputs]

        for ad, flat in zip(*gt.make_lists(adinputs, flat_list, force_ad=True)):
            # CJS: Since we're not saving the processed_arc before this, we
            # can't check for the tags. Instead, let's look for the WGT extn
            if not hasattr(ad[0], 'WGT'):
                log.warning("fitWavelength is only run on prepared GHOST arc"
                            " files - skipping {}".format(ad.filename))
                continue

            if self.timestamp_keys["extractProfile"] not in ad.phu.keywords:
                log.warning("extractProfile has not been run on {} - "
                            "skipping".format(ad.filename))
                continue

            if flat is None:
                log.warning("Could not find processed_flat calibration for "
                            "{} - skipping".format(ad.filename))
                continue

            try:
                poly_wave = self._get_polyfit_filename(ad, 'wavemod')
                poly_spat = self._get_polyfit_filename(ad, 'spatmod')
                poly_spec = self._get_polyfit_filename(ad, 'specmod')
                poly_rot = self._get_polyfit_filename(ad, 'rotmod')
                wpars = astrodata.open(poly_wave)
                spatpars = astrodata.open(poly_spat)
                specpars = astrodata.open(poly_spec)
                rotpars = astrodata.open(poly_rot)
            except IOError:
                log.warning("Cannot open required initial model files for {};"
                            " skipping".format(ad.filename))
                continue

            # CJS: line_list location is now in lookups/__init__.py
            arclinefile = os.path.join(os.path.dirname(polyfit_dict.__file__),
                                       line_list)
            arcwaves, arcfluxes = np.loadtxt(arclinefile, usecols=[1, 2]).T

            arm = GhostArm(arm=ad.arm(), mode=ad.res_mode())
            arm.spectral_format_with_matrix(flat[0].XMOD,
                                            wpars[0].data,
                                            spatpars[0].data,
                                            specpars[0].data,
                                            rotpars[0].data)

            extractor = Extractor(arm, None)  # slitview=None for this usage
            lines_out = extractor.find_lines(ad[0].data, arcwaves, inspect=False)

            #import pdb; pdb.set_trace()
            fitted_params, wave_and_resid = arm.read_lines_and_fit(
                wpars[0].data, lines_out)

            # CJS: Append the WFIT as an extension. It will inherit the
            # header from the science plane (including irrelevant/wrong
            # keywords like DATASEC) but that's not really a big deal.
            ad[0].WFIT = fitted_params

            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def flatCorrect(self, adinputs=None, **params):
        """
        Flat-correct an extracted GHOST object profile by extracting the
        profile from the relevant flat field using the object's extracted
        weights, and then perform simple division.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        flat: str/None
            Name of the (processed) standard flat to use for flat profile
            extraction. If None, the primitive will attempt to pull a flat
            from the calibrations database (or, if specified, the
            --user_cal processed_flat command-line option)
        slit: str/None
            Name of the (processed & stacked) slit image to use for extraction
            of the profile. If not provided/set to None, the primitive will
            attempt to pull a processed slit image from the calibrations
            database (or, if specified, the --user_cal processed_slit
            command-line option)
        slitflat: str/None
            Name of the (processed) slit flat image to use for extraction
            of the profile. If not provided, set to None, the RecipeSystem
            will attempt to pull a slit flat from the calibrations system (or,
            if specified, the --user_cal processed_slitflat command-line
            option)
        writeResult: bool
            Denotes whether or not to write out the result of profile
            extraction to disk. This is useful for both debugging, and data
            quality assurance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # CJS: See extractProfile() refactoring for explanation of changes
        slit_list = params["slit"]
        if slit_list is None:
            self.getProcessedSlit(adinputs)
            slit_list = [self._get_cal(ad, 'processed_slit')
                         for ad in adinputs]

        # CJS: I've renamed flat -> slitflat and obj_flat -> flat because
        # that's what the things are called! Sorry if I've overstepped.
        slitflat_list = params["slitflat"]
        if slitflat_list is None:
            self.getProcessedSlitFlat(adinputs)
            slitflat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        flat_list = params["flat"]
        if flat_list is None:
            self.getProcessedFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_flat')
                         for ad in adinputs]

        # TODO: Have gt.make_lists handle multiple auxiliary lists?
        _, slit_list = gt.make_lists(adinputs, slit_list, force_ad=True)
        _, slitflat_list = gt.make_lists(adinputs, slitflat_list, force_ad=True)
        _, flat_list = gt.make_lists(adinputs, flat_list, force_ad=True)

        for ad, slit, slitflat, flat, in zip(adinputs, slitflat_list,
                                        flat_list, flat_list):
            log.info(ad.info())

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by flatCorrect".
                            format(ad.filename))
                continue

            # CJS: failure to find a suitable auxiliary file (either because
            # there's no calibration, or it's missing) places a None in the
            # list, allowing a graceful continuation.
            if slit is None or slitflat is None or flat is None:
                log.warning("Unable to find calibrations for {}; "
                            "skipping".format(ad.filename))
                continue

            try:
                poly_wave = self._get_polyfit_filename(ad, 'wavemod')
                poly_spat = self._get_polyfit_filename(ad, 'spatmod')
                poly_spec = self._get_polyfit_filename(ad, 'specmod')
                poly_rot = self._get_polyfit_filename(ad, 'rotmod')
                wpars = astrodata.open(poly_wave)
                spatpars = astrodata.open(poly_spat)
                specpars = astrodata.open(poly_spec)
                rotpars = astrodata.open(poly_rot)
            except IOError:
                log.warning("Cannot open required initial model files for {};"
                            " skipping".format(ad.filename))
                continue

            res_mode = ad.res_mode()
            arm = GhostArm(arm=ad.arm(), mode=res_mode, 
                           detector_x_bin = ad.detector_x_bin(),
                           detector_y_bin = ad.detector_y_bin())
            arm.spectral_format_with_matrix(flat[0].XMOD,
                                            wpars[0].data,
                                            spatpars[0].data,
                                            specpars[0].data,
                                            rotpars[0].data)

            sview = SlitView(slit[0].data, slitflat[0].data, mode=res_mode)

            extractor = Extractor(arm, sview)
            extracted_flux, extracted_var = extractor.two_d_extract(
                arm.bin_data(flat[0].data), extraction_weights=ad[0].WGT)

            # Normalised extracted flat profile
            med = np.median(extracted_flux)
            extracted_flux /= med
            extracted_var /= med**2

            flatprof_ad = deepcopy(ad)
            flatprof_ad.update_filename(suffix='_extractedFlatProfile',
                                        strip=True)
            flatprof_ad[0].reset(extracted_flux, mask=None,
                                 variance=extracted_var)
            if params["write_result"]:
                flatprof_ad.write(overwrite=True)

            # Divide the flat field through the science data
            # Arithmetic propagates VAR correctly
            ad /= flatprof_ad

        return adinputs

    def rejectCosmicRays(self, adinputs=None, **params):
        """
        Reject cosmic rays from GHOST data.
        
        Parameters
        ----------
        n_steps: int
            The number of iterations that the LACosmic algorithm will make.
        subsampling: int
            The image subsampling factor LACosmic will use to generate the
            input images for the algorithm. There is really no reason to
            change this value from the default.
        sigma_lim: float
            The sigma-clipping limit to be applied to the noise map.
        f_lim: float
            The clipping limit for the fine-structure image.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        n_steps = params["n_steps"]
        subsampling = params["subsampling"]
        sigma_lim = params["sigma_lim"]
        f_lim = params["f_lim"]

        # Define the Laplacian and growth kernels for L.A.Cosmic
        laplace_kernel = np.array([
            [0.0, -1.0, 0.0],
            [-1.0, 4.0, -1.0],
            [0.0, -1.0, 0.0],
        ])
        growth_kernel = np.ones((3, 3), dtype=np.float64)

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by rejectCosmicRays".
                            format(ad.filename))
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

            log.stdinfo("Doing CR removal for {}".format(ad.filename))

            for ext in ad:
                # CJS: Added forced creation of DQ plane
                if ext.mask is None:
                    ext.mask = np.zeros_like(ext.data, dtype=np.uint16)
                log.stdinfo('-----')
                log.stdinfo("EXTVER {}".format(ext.hdr['EXTVER']))
                log.stdinfo('-----')

                # Define an array that will hold the cosmic ray flagging
                # Note that we're deliberately not using the BPM at this stage,
                # otherwise the algorithm will start searching for cosmic rays
                # around pixels that have been flagged bad for another reason.
                cosmic_bpm = np.zeros_like(ext.data, dtype=np.int16)

                # Start with a fresh copy of the data
                # Use numpy NaN to cover up any data detected bad so far
                # (i.e. 0 < BPM < 8)
                clean_data = np.copy(ext.data)
                clean_data[ext.mask > 0] = np.nan

                no_passes = 0
                new_crs = 1
                new_cr_pix = None

                while new_crs > 0 and no_passes < n_steps:
                    no_passes += 1
                    curr_crs = np.count_nonzero(cosmic_bpm)
                    if curr_crs > 0 and new_cr_pix is not None:
                        # Median out the pixels already defined as cosmic rays
                        log.stdinfo('Pass {}: Median over previously '
                                    'found CR pix'.format(no_passes))

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
                    log.stdinfo('Pass {}: Building sky model'.format(no_passes))
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
                    log.stdinfo('Pass {}: Computing Laplacian'.format(no_passes))
                    data_shape = ext.data.shape
                    # log.stdinfo(
                    #     'data array size: %s' % str(data_shape))
                    subsampl_data = np.repeat(np.repeat(
                        ext.data, subsampling, axis=1),
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
                    log.stdinfo('Pass {}: Constructing sigma map'.format(no_passes))
                    gain = ext.gain()
                    read_noise = ext.read_noise()
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
                    log.stdinfo('Pass {}: Flagging cosmic rays'.format(no_passes))
                    # Construct the fine-structure image
                    # (F, eqn 14 of van Dokkum)
                    m3 = scipy.ndimage.median_filter(subbed_data, size=[3, 3],
                                            mode='constant', cval=np.nan)
                    fine_struct = m3 - scipy.ndimage.median_filter(m3,
                                size=[7, 7], mode='constant', cval=np.nan)
                    # Pixels are flagged as being cosmic rays if:
                    # - The sig_detrend image (S') is > sigma_lim
                    # - The contrast between the Laplacian image (L+) and the
                    #   fine-structure image (F) is greater than f_lim
                    new_cr_pix = np.logical_and(sig_detrend > sigma_lim,
                                            (conv_data/fine_struct) > f_lim)
                    cosmic_bpm[new_cr_pix] = np.uint16(DQ.cosmic_ray)
                    new_crs = np.count_nonzero(cosmic_bpm) - curr_crs
                    log.stdinfo('Pass {}: Found {} CR pixels'.format(no_passes,
                                                                     new_crs))

                # For the moment, go with Mike Ireland's suggestion to require
                # a BPM update
                ext.mask |= cosmic_bpm

                log.debug('Flagged pix in BPM: {}'.format(
                    np.count_nonzero(ext.mask)))

            # CJS: Added this because you check for the keyword in this primitive!
            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
        return adinputs

    def responseCorrect(self, adinputs=None, **params):
        """
        Perform a sigma-clipping on the input data frame, such that any pixels
        outside the sigma threshold have their BPM value updated

        Parameters
        ----------
        skip: bool
            If True, this primitive will just return the adinputs immediately
        """
        if params.get('skip'):
            log.stdinfo('Skipping the response correct stage')
            return adinputs

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by responseCorrect".
                            format(ad.filename))
                continue

            # Timestamp; DO NOT update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def standardizeStructure(self, adinputs=None, **params):
        """
        The Gemini-level will try to attach an MDF because a GHOST image is
        tagged as SPECT. Rather than set parameters for that primitive to
        stop it from doing so, just override with a no-op primitive.
        
        Note: This could go in primitives_ghost.py if the SLITV version
        also no-ops.
        """
        return adinputs

    # CJS: Primitive has been renamed for consistency with other instruments
    # The geometry_conf.py file is not needed; all you're doing is tiling
    # extensions according to their DETSEC keywords, without gaps or rotations
    # so this shouldn't need any extra information.
    def tileArrays(self, adinputs=None, **params):
        """
        This primitive will tile the SCI frames of the input images, along
        with the VAR and DQ frames if they exist.
        """
        def simple_mosaic_function(ad):
            """
            This will go into MosaicAD as the default function.
            Being discussed within the SUSD team.
            """
            from gempy.mosaic.mosaicData import MosaicData
            from gempy.mosaic.mosaicGeometry import MosaicGeometry

            # Calling trim_to_data_section() corrects the WCS if the overscan
            # regions haven't been trimmed yet
            ad = gt.trim_to_data_section(ad, keyword_comments=self.keyword_comments)

            md = MosaicData()  # Creates an empty object
            md.data_list = []  # Not needed

            x_bin = ad.detector_x_bin()
            y_bin = ad.detector_y_bin()
            detsecs = [(k[0]//x_bin, k[1]//x_bin, k[2]//y_bin, k[3]//y_bin)
                       for k in ad.detector_section()]
            # One output block
            md.coords = {'amp_mosaic_coord': detsecs,
                         'amp_block_coord': detsecs}
            nxout = max(k[1] for k in detsecs)
            nyout = max(k[3] for k in detsecs)
            mg = MosaicGeometry({'blocksize': (nxout, nyout),
                                 'mosaic_grid': (1,1)})
            return md, mg

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        adoutputs = []
        for ad in adinputs:
            mo = MosaicAD(ad, mosaic_ad_function=simple_mosaic_function)
            ad_mos = mo.as_astrodata(tile=True)

            gt.mark_history(ad_mos, primname=self.myself(),
                            keyword=timestamp_key)
            ad_mos.update_filename(suffix=params["suffix"],
                                   strip=True)
            adoutputs.append(ad_mos)

            ad_mos.write(overwrite=True)
            # ad_mos.write(overwrite=True)

        return adoutputs

    # validateData() removed since inherited Standardize method will handle it

    def _get_polyfit_filename(self, ad, caltype):
        """
        Gets the filename of the relevant initial polyfit file for this
        input GHOST science image

        Returns
        -------
        str/None: Filename (including path) of the required polyfit file
        """
        log = self.log
        polyfit_dir = os.path.join(os.path.dirname(polyfit_dict.__file__),
                                   'Polyfit')

        # CJS: This is a method that only exists *if* the input is of type
        # GHOST, so no need to check
        arm = ad.arm()
        res_mode = ad.res_mode()
        key = 'GHOST_1_1_{}_{}'.format(arm, res_mode)

        try:
            poly_dict = getattr(polyfit_dict, '{}_dict'.format(caltype))
        except AttributeError:
            log.warning("Invalid polyfit calibration type ({}) requested for "
                        "{}".format(caltype, ad.filename))
            return None

        dates_avail = set([k.split('_')[-1] for k in poly_dict.keys()])
        # Safe to assume instrument won't be used after 2099...
        dates_avail = map(lambda x: datetime.strptime('20{}'.format(x),
                                            '%Y%m%d').date(), dates_avail)
        dates_avail.sort()

        # Determine the latest data that precedes the observing date
        date_obs = ad.ut_date()
        try:
            date_req = [_ for _ in dates_avail if _ <= date_obs][-1]
        except IndexError:
            log.warning("No polyfit available for {}".format(ad.filename))
            return None
        key += '_{}'.format(date_req.strftime('%y%m%d'))

        polyfit_file = poly_dict[key]
        # Prepend standard path if the filename doesn't start with '/'
        return polyfit_file if polyfit_file.startswith(os.path.sep) else \
            os.path.join(polyfit_dir, polyfit_file)

    def _compute_barycentric_correction(self, ad, return_wavl=True,
                                        loc=GEMINI_SOUTH_LOC):
        """
        Compute the baycentric correction factor for a given observation and
        location on the Earth.

        The barycentric correction compensates for (a) the motion of the Earth
        around the Sun, and (b) the motion of the Earth's surface due to
        planetary rotation. It can be returned as a line velocity correction,
        or a multiplicative factor with which to correct the wavelength scale;
        the default is the latter.

        The correction will be computed for all extensions of the input
        AstroData object.

        This method is built using astropy v2, and is developed from:
        https://github.com/janerigby/jrr/blob/master/barycen.py


        Parameters
        ----------
        ad : astrodata.AstroData
            The astrodata object from which to extract time and
            location information. If the ad is multi-extension, a correction
            factor will be returned for each extension.
        return_wavl : bool
            Denotes whether to return the correction as a wavelength
            correction factor (True) or a velocity (False). Defaults to True.

        Returns
        -------
        corr_facts: list of float
            The barycentric correction values, one per extension of the input
            ad.
        """

        # Set up a SkyCoord for this ad
        sc = astrocoord.SkyCoord(ad.phu.get('RA'), ad.phu.get('DEC'),
                                 unit=(u.deg, u.deg, ))

        # Compute central time of observation
        dt_start = datetime.combine(
            datetime.strptime(ad.phu.get('DATE-OBS'), '%Y-%m-%d').date(),
            datetime.strptime(ad.phu.get('UTSTART'), '%H:%M:%S.%f').time(),
        )

        corr_facts = []
        for ext in ad:
            dt_midp = dt_start + timedelta(
                seconds=ext.hdr.get('EXPTIME')/2.0
            )
            dt_midp = Time(dt_midp)
            # ICRS position & vel of Earth geocenter
            ep, ev = astrocoord.solar_system.get_body_barycentric_posvel(
                'earth', dt_midp
            )
            # GCRS position & vel of observatory (loc)
            op, ov = loc.get_gcrs_posvel(dt_midp)
            # Velocities can be simply added (are axes-aligned)
            vel = ev + ov

            # Get unit ICRS vector in direction of observation
            sc_cart = sc.icrs.represent_as(
                astrocoord.UnitSphericalRepresentation
            ).represent_as(
                astrocoord.CartesianRepresentation
            )

            corr_fact = sc_cart.dot(vel).to(u.km/u.s)
            if return_wavl:
                corr_fact = 1.0 + (corr_fact / const.c)

            corr_facts.append(corr_fact)

        return corr_facts
