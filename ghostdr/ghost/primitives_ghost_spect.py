#
#                                                                  gemini_python
#
#                                                      primitives_ghost_spect.py
# ------------------------------------------------------------------------------
import os
import numpy as np
import math
from copy import deepcopy
import scipy
import scipy.signal as signal
from scipy.optimize import leastsq
import functools
from datetime import datetime, date, time, timedelta
import re
import astropy.coordinates as astrocoord
import astropy.io.fits as astropyio
from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astropy.stats import sigma_clip
from scipy import interpolate
from pysynphot import observation, spectrum

import astrodata

from geminidr.gemini.lookups import DQ_definitions as DQ

from gempy.gemini import gemini_tools as gt
from gempy.mosaic.mosaicAD import MosaicAD

from .polyfit import GhostArm, Extractor, SlitView
from .polyfit.ghost import GhostArm

from .primitives_ghost import GHOST, filename_updater
from .parameters_ghost_spect import ParametersGHOSTSpect

from .lookups import polyfit_dict, line_list, keyword_comments, targetn_dict

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

GEMINI_SOUTH_LOC = astrocoord.EarthLocation.from_geodetic((-70, 44, 12.096),
                                                          (-30, 14, 26.700),
                                                          height=2722.,
                                                          ellipsoid='WGS84')


@parameter_override
class GHOSTSpect(GHOST):
    """
    Primitive class for processing GHOST science data.

    This class contains the primitives necessary for processing GHOST science
    data, as well as all related calibration files from the main spectrograph
    cameras. Slit viewer images are processed with another primitive class
    (:class:`ghostdr.ghost.primitives_ghost_slit.GHOSTSlit`).
    """

    """Applicable tagset"""
    tagset = set(["GEMINI", "GHOST"])  # NOT SPECT because of bias/dark

    def __init__(self, adinputs, **kwargs):
        super(GHOSTSpect, self).__init__(adinputs, **kwargs)
        self.parameters = ParametersGHOSTSpect
        self.keyword_comments.update(keyword_comments.keyword_comments)

    def addWavelengthSolution(self, adinputs=None, **params):
        """
        Compute and append a wavelength solution for the data.

        The GHOST instrument is designed to be very stable over a long period
        of time, so it is not strictly necessary to take arcs for every
        observation. The alternative is use the arcs taken most recently
        before and after the observation of interest, and compute an
        average of their wavelength solutions.

        The average is weighted by
        the inverse of the time between each arc observation and science
        observation. E.g., if the 'before' arc is taken 12 days before the
        science observation, and the 'after' arc is taken 3 days after the
        science observation, then the 'after' arc will have a weight of 80%
        in the final wavelength solution (12/15), and the 'before' arc 20%
        (3/15).

        In the event that either a 'before' arc can't be found but an 'after'
        arc can, or vice versa, the wavelength solution from the arc that was
        found will be applied as-is. If neither a 'before' nor 'after' arc can
        be found, an IOError will be raised.

        It is possible to explicitly pass which arc files to use as
        the ``arc`` parameter. This should be a list of two-tuples, with each
        tuple being of the form
        ``('before_arc_filepath', 'after_arc_filepath')``. This list must be
        the same length as the list of ``adinputs``, with a one-to-one
        correspondence between the two lists.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        arc: list of two-tuples
            A list of two-tuples, with each tuple corresponding to an element of
            the ``adinputs`` list. Within each tuple, the two elements are the
            designated 'before' and 'after' arc for that observation.
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

            # self.getProcessedArc(ad, howmany=2)
            # if not found_arcs:
            #     try:
            #         arcs_calib = self._get_cal(ad, 'processed_arc', )
            #         log.stdinfo('Found following arcs: {}'.format(
            #             ', '.join([_ for _ in arcs_calib])
            #         ))
            #         arc_before, arc_after = self._get_cal(ad, 'processed_arc',)
            #     except (TypeError, ValueError):
            #         # Triggers if only one arc, or more than two
            #         arc_before = self._get_cal(ad, 'processed_arc',)[0]
            #         arc_after = None

            if not found_arcs:
                # Fetch the arc_before and arc_after in sequence
                arc_before = self._request_bracket_arc(ad, before=True)
                arc_after = self._request_bracket_arc(ad, before=False)

            if arc_before is None and arc_after is None:
                raise IOError('No valid arcs found for {}'.format(ad.filename))

            log.stdinfo('Arcs for {}: \n'
                        '   before: {}\n'
                        '    after: {}'.format(ad.filename,
                                               arc_before, arc_after))

            # Stand up a GhostArm instance for this ad
            gs = GhostArm(arm=ad.arm(), mode=ad.res_mode(),
                          detector_x_bin=ad.detector_x_bin(),
                          detector_y_bin=ad.detector_y_bin())

            if arc_before is None:
                # arc = arc_after
                arc_after = astrodata.open(arc_after)
                wfit = gs.evaluate_poly(arc_after[0].WFIT)
                ad.phu.set('ARCIM_A', os.path.abspath(arc_after.path),
                           "'After' arc image")
            elif arc_after is None:
                # arc = arc_before
                arc_before = astrodata.open(arc_before)
                wfit = gs.evaluate_poly(arc_before[0].WFIT)
                ad.phu.set('ARCIM_B', os.path.abspath(arc_before.path),
                           "'Before' arc image")
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
                ad.phu.set('ARCIM_A', os.path.abspath(arc_after.path),
                           self.keyword_comments['ARCIM_A'])
                ad.phu.set('ARCIM_B', os.path.abspath(arc_before.path),
                           self.keyword_comments['ARCIM_B'])
                ad.phu.set('ARCWT_A', weight_a,
                           self.keyword_comments['ARCWT_A'])
                ad.phu.set('ARCWT_B', weight_b,
                           self.keyword_comments['ARCWT_B'])

            # rebin the wavelength fit to match the rest of the extensions
            for _ in range(int(math.log(ad.detector_x_bin(), 2))):
                wfit = wfit[:, ::2] + wfit[:, 1::2]
                wfit /= 2.0

            ad[0].WAVL = wfit

            # FIXME Wavelength unit needs to be in output ad

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def applyFlatBPM(self, adinputs=None, **params):
        """
        Find the flat relevant to the file(s) being processed, and merge the
        flat's BPM into the target file's.

        GHOST does not use flat subtraction in the traditional sense; instead,
        the extracted flat profile is subtracted from the extracted object
        profile. This means that the BPM from the flat needs to be applied to
        the object file before profile extraction, and hence well before actual
        flat correction is performed.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        flat: str/None
            Name (full path) of the flatfield to use. If None, try:
        flatstream: str/None
            Name of the stream containing the flatfield as the first
            item in the stream. If None, the calibration service is used
        write_result: bool
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
            ) or flat.detector_y_bin() != ad.detector_y_bin():
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

            ad.phu.set('FLATBPM', os.path.abspath(flat.path),
                       self.keyword_comments['FLATBPM'])

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
            if params["write_result"]:
                ad.phu.set('PROCIMG', os.path.abspath(ad.path),
                           keyword_comments.keyword_comments['PROCIMG'])
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

            # Get or compute the correction factor
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
        Perform a sigma-clipping on the input data frame.

        Parameters
        ----------
        sigma: float/None
            The sigma value to be used for clipping.
        bpm_value: int/None
            The integer value to be applied to the data BPM where the sigma
            threshold is exceeded. Defaults to 1 (which is the generic bad
            pixel flag). Note that the final output BPM is made using a
            bitwise_or operation.
        iters : int/None
            Number of sigma clipping iterations to perform. Default is None,
            which will continue sigma clipping until no further points are
            masked.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        sigma = params["sigma"]
        bpm_value = params["bpm_value"]
        iters = params["iters"]

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by clipSigmaBPM".
                            format(ad.filename))
                continue

            for ext in ad:
                extver = ext.hdr['EXTVER']
                if ext.mask is not None:
                    # Perform the sigma clip
                    clipd = sigma_clip(ext.data, sigma=sigma,
                                       iters=iters, copy=True)
                    # Convert the mask from the return into 0s and 1s and
                    # bitwise OR into the ext BPM
                    clipd_mask = clipd.mask.astype(ext.mask.dtype)
                    ext.mask |= clipd_mask * bpm_value

                    log.stdinfo('   {}:{}: nPixMasked: {:9d} / {:9d}'.format(
                        ad.filename, extver, np.sum(clipd_mask), ext.data.size))

                    # Original implementaion
                    # mean_data = np.mean(ext.data)
                    # sigma_data = np.std(ext.data)
                    # mask_map = (np.abs(ext.data-mean_data) > sigma*sigma_data)
                    # if bpm_value:  # might call with None for diagnosis
                    #     ext.mask[mask_map] |= bpm_value
                    #
                    # log.stdinfo('   {}:{}: nPixMasked: {:9d} / {:9d}'.format(
                    #     ad.filename, extver, np.sum(mask_map), ext.data.size))
                else:
                    log.warning('No DQ plane in {}:{}'.format(ad.filename,
                                                              extver))

            # Timestamp; DO NOT update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        return adinputs

    def darkCorrect(self, adinputs=None, **params):
        """
        Dark-correct GHOST observations.

        This primitive, at its core, simply copies the standard
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

        adinputs_orig = list(adinputs)
        if isinstance(params.get('dark', None), list):
            params['dark'] = [params['dark'][i] for i in range(len(adinputs))
                              if not adinputs[i].phu.get(timestamp_key)]
        adinputs = [_ for _ in adinputs if not _.phu.get(timestamp_key)]
        if len(adinputs) != len(adinputs_orig):
            log.stdinfo('The following files have already been processed by '
                        'darkCorrect and will not be further modified: '
                        '{}'.format(', '.join([_.filename for _ in adinputs_orig
                                               if _ not in adinputs])))

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
            ad.phu.set('DARKIM', os.path.abspath(dark.path),
                       self.keyword_comments["DARKIM"])
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs_orig

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
        sky_correct: bool
            Denotes whether or not to correct for the sky profile during the
            object extraction. Defaults to True, although it should be altered
            to False when processing flats or arcs.
        writeResult: bool
            Denotes whether or not to write out the result of profile
            extraction to disk. This is useful for both debugging, and data
            quality assurance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # This primitive modifies the input AD structure, so it must now
        # check if the primitive has already been applied. If so, it must be
        # skipped.
        adinputs_orig = list(adinputs)
        adinputs = [_ for _ in adinputs if not _.phu.get(timestamp_key)]
        if len(adinputs) != len(adinputs_orig):
            log.stdinfo('extractProfile is skipping the following files, which '
                        'already have extracted profiles: '
                        '{}'.format(','.join([_.filename for _ in adinputs_orig
                                              if _ not in adinputs])))

        # CJS: Heavily edited because of the new AD way
        # Get processed slits, slitFlats, and flats (for xmod)
        # slits and slitFlats may be provided as parameters
        slit_list = params["slit"]
        # log.stdinfo('slit_list before processing:')
        # log.stdinfo('   {}'.format(slit_list))
        if slit_list is not None and isinstance(slit_list, list):
            slit_list = [slit_list[i] for i in range(len(slit_list))
                         if adinputs_orig[i] in adinputs]
        if slit_list is None:
            # CJS: This populates the calibrations cache (dictionary) with
            # "processed_slit" filenames for each input AD
            self.getProcessedSlit(adinputs)
            # This then gets those filenames
            slit_list = [self._get_cal(ad, 'processed_slit')
                         for ad in adinputs]
        # log.stdinfo('slit_list after processing:')
        # log.stdinfo('   {}'.format(slit_list))

        slitflat_list = params["slitflat"]
        if slitflat_list is not None and isinstance(slitflat_list, list):
            slitflat_list = [slitflat_list[i] for i in range(len(slitflat_list))
                             if adinputs_orig[i] in adinputs]
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
            # CJS: failure to find a suitable auxiliary file (either because
            # there's no calibration, or it's missing) places a None in the
            # list, allowing a graceful continuation.
            if slit is None or slitflat is None or flat is None:
                log.warning("Unable to find calibrations for {}; "
                            "skipping".format(ad.filename))
                continue

            # CJS: Changed to log.debug() and changed the output
            log.stdinfo("Slit parameters: ")
            log.stdinfo("   processed_slit: {}".format(slit.filename))
            log.stdinfo("   processed_slitflat: {}".format(slitflat.filename))
            log.stdinfo("   processed_flat: {}".format(flat.filename))

            res_mode = ad.res_mode()
            arm = GhostArm(arm = ad.arm(), mode = res_mode,
                           detector_x_bin = ad.detector_x_bin(),
                           detector_y_bin = ad.detector_y_bin())

            # CJS: Heavy refactor. Return the filename for each calibration
            # type. Eliminates requirement that everything be updated
            # simultaneously.
            # key = self._get_polyfit_key(ad)
            # log.stdinfo("Polyfit key selected: {}".format(key))
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
            _, _, extracted_weights = extractor.one_d_extract(
                ad[0].data, correct_for_sky=params['sky_correct'],
            )
            extracted_flux, extracted_var = extractor.two_d_extract(
                ad[0].data,
                extraction_weights=extracted_weights,
            )

            # CJS: Since you don't use the input AD any more, I'm going to
            # modify it in place, in line with your comment that you're
            # considering this.
            ad[0].reset(extracted_flux, mask=None, variance=extracted_var)
            ad[0].WGT = extracted_weights

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
            ad[0].hdr['DATADESC'] = ('Order-by-order processed science data',
                                     self.keyword_comments['DATADESC'])
            if params["write_result"]:
                ad.write(overwrite=True)

        return adinputs_orig

    def interpolateAndCombine(self, adinputs=None, **params):
        """
        Combine the independent orders from the input ADs into a single,
        over-sampled spectrum

        Parameters
        ----------
        scale : str
            Denotes what scale to generate for the final spectrum. Currently
            available are:
            ``'loglinear'``
            Default is ``'loglinear'``.
        skip : bool
            Set to ``True`` to skip this primitive. Defaults to ``False``.
        oversample : int or float
            The factor by which to (approximately) oversample the final output
            spectrum, as compared to the input spectral orders. Defaults to 2.
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

            if params['skip']:
                log.warning('Skipping interpolateAndCombine for {}'.format(
                    ad.filename
                ))
                continue

            # Determine the wavelength bounds of the file
            min_wavl, max_wavl = np.min(ad[0].WAVL), np.max(ad[0].WAVL)
            logspacing = np.median(
                np.log(ad[0].WAVL[:, 1:]) - np.log(ad[0].WAVL[:, :-1])
            )
            # Form a new wavelength scale based on these extremes
            if params['scale'] == 'loglinear':
                wavl_grid = np.exp(
                    np.linspace(np.log(min_wavl), np.log(max_wavl),
                                num=int(
                                    (np.log(max_wavl) - np.log(min_wavl)) /
                                    (logspacing / float(params['oversample']))
                                ))
                )
            else:
                raise ValueError('interpolateAndCombine does not understand '
                                 'the scale {}'.format(params['scale']))

            # Create a final spectrum and (inverse) variance to match
            # (One plane per object)
            spec_final = np.zeros(wavl_grid.shape + (3, ))
            var_final = np.inf * np.ones(wavl_grid.shape + (3, ))

            # Loop over each input order, making the output spectrum the
            # result of the weighted average of itself and the order
            # spectrum
            for order in range(ad[0].data.shape[0]):
                for ob in range(ad[0].data.shape[-1]):
                    log.stdinfo('Re-gridding order {:2d}, obj {:1d}'.format(
                        order, ob,
                    ))
                    flux_for_adding = np.interp(wavl_grid,
                                                ad[0].WAVL[order],
                                                ad[0].data[order, :, ob],
                                                left=0, right=0)
                    ivar_for_adding = np.interp(wavl_grid,
                                                ad[0].WAVL[order],
                                                1.0 /
                                                ad[0].variance[order, :, ob],
                                                left=0, right=0)
                    spec_comp, ivar_comp = np.ma.average(
                        np.asarray([spec_final[:, ob], flux_for_adding]),
                        weights=np.asarray([1.0 / var_final[:, ob],
                                            ivar_for_adding]),
                        returned=True, axis=0,
                    )
                    spec_final[:, ob] = deepcopy(spec_comp)
                    var_final[:, ob] = deepcopy(1.0 / ivar_comp)

            # import pdb;
            # pdb.set_trace()

            # MCW, 180501 - Keep initial data, append interp'd data

            ad_interp = deepcopy(ad)
            # Can't use .reset without looping through extensions
            ad_interp[0].data = spec_final
            ad_interp[0].variance = var_final
            ad_interp[0].WAVL = wavl_grid
            try:
                del ad_interp[0].WGT
            except AttributeError:
                pass
            ad_interp[0].hdr['DATADESC'] = ('Interpolated data',
                                            self.keyword_comments['DATADESC'], )
            ad.append(ad_interp[0])

            # Timestamp & update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def findApertures(self, adinputs=None, **params):
        """
        Locate the slit aperture, parametrized by a :any:`polyfit` model.

        The primitive locates the slit apertures within a GHOST frame,
        and inserts a :any:`polyfit` model into a new extension on each data
        frame. This model is placed into a new ``.XMOD`` attribute on the
        extension.
        
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
            condit = np.where(np.abs(
                (medfilt - data)/(medfilt+1)) > 200
                              ) and np.where(data > np.average(data))
            data[condit] = medfilt[condit]

            # Convolve the flat field with the slit profile
            flat_conv = ghost_arm.slit_flat_convolve(
                data,
                slit_profile=slitview.slit_profile(arm=arm),
                spatpars=spatpars[0].data, microns_pix=slitview.microns_pix,
                xpars=xpars[0].data
            )

            flat_conv = signal.medfilt2d(flat_conv, (5, 5))

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
        Fit wavelength solution to a GHOST ARC frame.

        This primitive should only be applied to a reduce GHOST ARC frame. Any
        other files passed through this primitive will be skipped.

        The primitive will use the arc line files stored in the same location
        as the initial :any:`polyfit` models kept in the ``lookups`` system.

        This primitive uses no special parameters.
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

            if self.timestamp_keys["extractProfile"] not in ad.phu:
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
            # import pdb;pdb.set_trace()
            lines_out = extractor.find_lines(ad[0].data, arcwaves,
                                             inspect=False)
            
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
        Flat-correct an extracted GHOST profile using a flat profile.

        This primitive works by extracting the
        profile from the relevant flat field using the object's extracted
        weights, and then performs simple division.

        .. warning::
            While the primitive is working, it has been found that the
            underlying algorithm is flawed. A new algorithm is being developed.

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
        sfx = params["suffix"]

        adinputs_orig = list(adinputs)
        adinputs = [_ for _ in adinputs if not _.phu.get(timestamp_key)]
        if len(adinputs) != len(adinputs_orig):
            log.stdinfo('flatCorrect is skipping the following files, '
                        'which are already flat corrected: '
                        '{}'.format(','.join([_ for _ in adinputs_orig
                                              if _ not in adinputs])))

        # CJS: See extractProfile() refactoring for explanation of changes
        slit_list = params["slit"]
        if slit_list is not None and isinstance(slit_list, list):
            slit_list = [slit_list[i] for i in range(len(slit_list))
                         if adinputs_orig[i] in adinputs]
        if slit_list is None:
            self.getProcessedSlit(adinputs)
            slit_list = [self._get_cal(ad, 'processed_slit')
                         for ad in adinputs]

        # CJS: I've renamed flat -> slitflat and obj_flat -> flat because
        # that's what the things are called! Sorry if I've overstepped.
        slitflat_list = params["slitflat"]
        if slitflat_list is not None and isinstance(slitflat_list, list):
            slitflat_list = [slitflat_list[i] for i in range(len(slitflat_list))
                         if adinputs_orig[i] in adinputs]
        if slitflat_list is None:
            self.getProcessedSlitFlat(adinputs)
            slitflat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        flat_list = params["flat"]
        if flat_list is not None and isinstance(flat_list, list):
            flat_list = [flat_list[i] for i in range(len(flat_list))
                         if adinputs_orig[i] in adinputs]
        if flat_list is None:
            self.getProcessedFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_flat')
                         for ad in adinputs]

        # TODO: Have gt.make_lists handle multiple auxiliary lists?
        _, slit_list = gt.make_lists(adinputs, slit_list, force_ad=True)
        _, slitflat_list = gt.make_lists(adinputs, slitflat_list, force_ad=True)
        _, flat_list = gt.make_lists(adinputs, flat_list, force_ad=True)

        for ad, slit, slitflat, flat, in zip(adinputs, slit_list,
                                             slitflat_list, flat_list):
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
                           detector_x_bin= ad.detector_x_bin(),
                           detector_y_bin= ad.detector_y_bin()
                           )
            arm.spectral_format_with_matrix(flat[0].XMOD,
                                            wpars[0].data,
                                            spatpars[0].data,
                                            specpars[0].data,
                                            rotpars[0].data,
                                            )

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
                # Record this as the flat profile used
                ad.phu.set('FLATPROF', os.path.abspath(flatprof_ad.path),
                           self.keyword_comments['FLATPROF'])
            ad.phu.set('FLATIMG', os.path.abspath(flat.path),
                       keyword_comments.keyword_comments['FLATIMG'])
            ad.phu.set('SLITIMG', os.path.abspath(slit.path),
                       keyword_comments.keyword_comments['SLITIMG'])
            ad.phu.set('SLITFLAT', os.path.abspath(slitflat.path),
                       keyword_comments.keyword_comments['SLITFLAT'])

            # Divide the flat field through the science data
            # Arithmetic propagates VAR correctly
            ad /= flatprof_ad

            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs_orig

    def formatOutput(self, adinputs=None, **params):
        """
        Generate an output FITS file containing the data requested by the user.

        This primitive should not be called until *all* required
        processing steps have been performed on the data. THe resulting FITS
        file cannot be safely passed through to other primitives.

        .. note::
            All of the extra data packaged up by this primitive can also be
            obtained by using the ``write_result=True`` flag on selected
            other primitives. ``formatOutput`` goes and finds those output
            files, and then packages them into the main output file for
            convenience.

        Parameters
        ----------
        detail: str
            The level of detail the user would like in their final output file.

            Note that, in order to preserve the ordering of FITS file
            extensions, the options are sequential; each option will
            provide all the data of less-verbose options.

            Valid options are:

            ``default``
                Only returns the extracted, fully-processed object(s) and sky
                spectra. In effect, this causes ``formatOutput`` to do nothing.
                This includes computed variance data for each plane.

            ``processed_image``
                The option returns the data that have been bias and dark
                corrected, and has the flat BPM applied (i.e. the state the
                data are in immediately prior to profile extraction).

            ``flat_profile``
                This options includes the extracted flat profile used for
                flat-fielding the data.

            ``sensitivity_curve``
                This option includes the sensitivity calculated at the
                ``responseCorrect`` step of reduction.
        """

        # This should be the list of allowed detail descriptors in order of
        # increasing verbosity
        ALLOWED_DETAILS = ['default', 'processed_image', 'flat_profile',
                           'sensitivity_curve', ]

        log = self.log

        timestamp_key = self.timestamp_keys[self.myself()]
        sfx = params['suffix']

        if params['detail'] not in ALLOWED_DETAILS:
            raise ValueError('formatOutput: detail option {} not known. '
                             'Please use one of: {}'.format(
                params['detail'],
                ', '.join(ALLOWED_DETAILS),
            ))

        detail_index = ALLOWED_DETAILS.index(params['detail'])

        for ad in adinputs:
            # Move sequentially through the various levels of detail, adding
            # them as we go along
            # ad[0].hdr['DATADESC'] = ('Fully-reduced data',
            #                          self.keyword_comments['DATADESC'], )

            if ALLOWED_DETAILS.index('processed_image') <= detail_index:
                # Locate the processed image data
                fn = ad.phu.get('PROCIMG', None)
                if fn is None:
                    raise RuntimeError('The processed image file name for {} '
                                       'has not been '
                                       'recorded'.format(ad.filename))
                try:
                    proc_image = astrodata.open(fn)
                except astrodata.AstroDataError:
                    raise RuntimeError('You appear not to have written out '
                                       'the result of image processing to '
                                       'disk.')

                log.stdinfo('Opened processed image file {}'.format(fn))
                ad.append(proc_image[0])
                ad[-1].hdr['DATADESC'] = ('Processed image',
                                          self.keyword_comments['DATADESC'])

            if ALLOWED_DETAILS.index('flat_profile') <= detail_index:
                # Locate the flat profile data
                fn = ad.phu.get('FLATPROF', None)
                if fn is None:
                    raise RuntimeError('The flat profile file name for {} '
                                       'has not been '
                                       'recorded'.format(ad.filename))
                try:
                    proc_image = astrodata.open(fn)
                except astrodata.AstroDataError:
                    raise RuntimeError('You appear not to have written out '
                                       'the result of flat profiling to '
                                       'disk.')

                log.stdinfo('Opened flat profile file {}'.format(fn))
                # proc_image[0].WGT = None
                try:
                    del proc_image[0].WGT
                except AttributeError:
                    pass
                ad.append(proc_image[0])
                ad[-1].hdr['DATADESC'] = ('Flat profile',
                                          self.keyword_comments['DATADESC'])

            if ALLOWED_DETAILS.index('sensitivity_curve') <= detail_index:
                fn = ad.phu.get('SENSFUNC', None)
                if fn is None:
                    raise RuntimeError('The sensitivity curve file name for {} '
                                       'has not been '
                                       'recorded'.format(ad.filename))
                try:
                    proc_image = astrodata.open(fn)
                except astrodata.AstroDataError:
                    raise RuntimeError('You appear not to have written out '
                                       'the result of sensitivity calcs to '
                                       'disk.')

                log.stdinfo('Opened sensitivity curve file {}'.format(fn))
                # proc_image[0].WGT = None
                try:
                    del proc_image[0].WGT
                except AttributeError:
                    pass
                try:
                    del proc_image[0].WAVL
                except AttributeError:
                    pass
                ad.append(proc_image[0])
                ad[-1].hdr['DATADESC'] = ('Sensitivity curve (blaze func.)',
                                          self.keyword_comments['DATADESC'])

            # import pdb; pdb.set_trace();

            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=sfx, strip=True)
            ad.write(overwrite=True)

        return adinputs

    def rejectCosmicRays(self, adinputs=None, **params):
        """
        Reject cosmic rays from GHOST data.

        .. warning::
            This primitive is now deprecated - cosmic ray rejection is now
            handled as part of the profile extraction process.
        
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
        raise DeprecationWarning('Cosmic ray rejections is now handled '
                                 'as part of the profile extraction process.')

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
                        # log.stdinfo('Padded array size: %s' %
                        #             str(pad_data.shape))
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
                    log.stdinfo('Pass {}: Computing Laplacian'.format(
                        no_passes)
                    )
                    data_shape = ext.data.shape
                    # log.stdinfo(
                    #     'data array size: %s' % str(data_shape))
                    subsampl_data = np.repeat(np.repeat(
                        ext.data, subsampling, axis=1),
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
                    log.stdinfo('Pass {}: Constructing sigma map'.format(
                        no_passes
                    ))
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
                    log.stdinfo('Pass {}: Flagging cosmic rays'.format(
                        no_passes
                    ))
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

            # CJS: Added this because you check for the keyword in
            # this primitive!
            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def responseCorrect(self, adinputs=None, **params):
        """
        Use a standard star observation and reference spectrum to provide
        absolute flux calibration.

        This primitive follows the basic pattern for determining absolute flux
        from an observed standard with a relative flux scale (e.g. counts) and
        an absolute flux-calibrated reference spectrum:

        - Dividing the standard star observation (in counts or electrons per
          pixel) by
          the exposure time (in s), and then by the standard star reference
          spectrum (in some unit of flux, e.g. erg/cm^2/s/A) gives a
          sensitivity curve in units of, in this example, counts / erg.
        - Dividing the object spectrum by the exposure time (i.e. converting 
          to counts per pixel per second) by the sensitivity curve
          (counts / flux unit) yields the object spectrum in the original flux
          units of the standard star reference spectrum.

        Parameters
        ----------
        skip: bool
            If True, this primitive will just return the adinputs immediately
        std : str, giving a relative or absolute file path
            The name of the reduced standard star observation. Defaults to
            None, at which point a ValueError is thrown.
        std_spec: str, giving a relative or absolute file path
            The name of the file where the standard star spectrum (the
            reference, not the observed one) is stored. Defaults to None,
            at which point a fatal error will be thrown.

            Spectral standard references should be in the format provided
            by Space Telescope Science Institute, e.g., from
            ftp://ftp.stsci.edu/cdbs/current_calspec/. If the standard reference
            is taken from elsewhere, it needs to obey the following
            format rules:

            - The reference data is in the first science extension of the FITS
              file;
            - The data must be in FITS table format, with columns named
              ``'FLUX'`` and ``'WAVELENGTH'``;
            - The first science extension must have a header card named
              ``'TUNIT2'``, which should contain the FITS-compatible
              flux unit name corresponding to the data in the ``'FLUX'`` column.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        if params.get('skip'):
            log.stdinfo('Skipping the response (standard star) correction '
                        'step')
            return adinputs

        if params['std'] is None:
            raise ValueError('No standard star provided')

        # Let the astrodata routine handle any issues with actually opening
        # the FITS file
        std = astrodata.open(params['std'])

        # Need to find the reference standard star spectrum
        # Use the one passed by the user in the first instance, otherwise
        # attempt to locate a remote one
        # Throw an error if none found
        if params['std_spec']:
            # TODO Will need to be modified to use Gemini service
            std_spec = astropyio.open(params['std_spec'])
            bunit = std_spec[1].header['TUNIT2']
        else:
            raise ValueError('No standard reference spectrum found/supplied')

        # Re-grid the standard reference spectrum onto the wavelength grid of
        # the observed standard
        regrid_std_ref = np.zeros(std[0].data.shape[:-1], dtype=np.float32)
        for od in range(std[0].data.shape[0]):
            regrid_std_ref[od] = self._regrid_spect(
                std_spec[1].data['FLUX'],
                std_spec[1].data['WAVELENGTH'],
                std[0].WAVL[od, :],
                waveunits='angstrom'
            )

        # Figure out which object is actually the standard observation
        # (i.e. of the dimensions [order, wavl, object], figure which of the
        # three objects is actually the spectrum (another will be sky, and
        # the third probably empty)
        objn = targetn_dict.targetn_dict['object']
        target = -1
        if std.phu['TARGET1'] == objn: target = 0
        if std.phu['TARGET2'] == objn: target = 1
        if target < 0:
            raise ValueError(
                'Cannot determine which IFU contains standard star spectrum.'
            )

        # Compute the sensitivity function
        sens_func = (std[0].data[:, :, target] /
                     std[0].hdr['EXPTIME']) / regrid_std_ref
        sens_func_var = (std[0].variance[:, :, target] /
                         std[0].hdr['EXPTIME']**2) / regrid_std_ref**2

        # MCW 180501
        # The sensitivity function requires significant smoothing in order to
        # prevent noise from the standard being transmitted into the data
        # The easiest option is to perform a parabolic curve fit to each order
        # QUADRATIC
        # fitfunc = lambda p, x: p[0] + p[2] * ((x - p[1])**2)
        # LINEAR
        fitfunc = lambda p, x: p[0] + (p[1] * x)
        errfunc = lambda p, x, y, yerr: np.abs(fitfunc(p, x) - y) / np.sqrt(yerr)
        # import pdb; pdb.set_trace();
        sens_func_fits = [
            p for p, success in [leastsq(errfunc,
                                         # QUADRATIC
                                         # [np.average(sens_func[od],
                                         #             weights=1./np.sqrt(
                                         #                 sens_func_var[od])),
                                         #  np.median(std[0].WAVL[od, :]),
                                         #  1.0],
                                         # LINEAR
                                         [np.average(sens_func[od, :],
                                          weights=1. / np.sqrt(
                                              sens_func_var[od])),
                                          0.],
                                         args=(std[0].WAVL[od, :],
                                               sens_func[od, :],
                                               sens_func_var[od, :])
                                         )
                                 for od in range(sens_func.shape[0])
                                 ]
            # if success
        ]

        # import pdb; pdb.set_trace();

        for ad in adinputs:

            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by responseCorrect".
                            format(ad.filename))
                continue

            # Check that the ad matches the standard
            if ad.res_mode() != std.res_mode():
                raise ValueError('Resolution modes do not match for '
                                 '{} and {}'.format(ad.filename, std.filename))
            if ad.arm() != std.arm():
                raise ValueError('Spectrograph arms do not match for '
                                 '{} and {}'.format(ad.filename, std.filename))
            if ad.detector_y_bin() != std.detector_y_bin() or \
                    ad.detector_x_bin() != std.detector_x_bin():
                raise ValueError('Binning does not match for '
                                 '{} and {}'.format(ad.filename, std.filename))

            # Interpolate the sensitivity function onto the wavelength
            # grid of this ad
            # Note that we can get away with this instead of a more
            # in-depth, flux-conserving regrid because:
            # (a) The sensitivity curve units do not depend on wavelength;
            # (b) The wavelength shifts involved are very small
            sens_func_regrid = np.zeros(ad[0].data.shape, dtype=np.float32)
            # sens_func_regrid_var = np.inf * np.ones(ad[0].data.shape,
            #                                         dtype=np.float32)
            for ob in range(ad[0].data.shape[-1]):
                for od in range(ad[0].data.shape[0]):
                    # import pdb; pdb.set_trace();
                    sens_func_regrid[od, :, ob] = fitfunc(
                        sens_func_fits[od], ad[0].WAVL[od, :]
                    )
                    # if od == 29:
                    #     import pdb; pdb.set_trace();
                    # sens_func_regrid[od, :, ob] = np.interp(
                    #     ad[0].WAVL[od, :],
                    #     std[0].WAVL[od, :],
                    #     sens_func[od, :],
                    #     left=0, right=0,
                    # )
                    # sens_func_regrid_var[od, :, ob] = np.interp(
                    #     ad[0].WAVL[od, :],
                    #     std[0].WAVL[od, :],
                    #     sens_func_var[od, :],
                    #     left=0, right=0,
                    # )


            # Easiest way to response correct is to stand up a new AstroData
            # instance containing the sensitivity function - this will
            # automatically handle, e.g., the VAR re-calculation
            sens_func_ad = deepcopy(ad)
            sens_func_ad.update_filename(suffix='_sensFunc', strip=True)
            sens_func_ad[0].data = sens_func_regrid
            # sens_func_ad[0].variance = sens_func_regrid_var
            # sens_func_ad[0].variance = None
            try:
                sens_func_ad[0].variance = None
            except AttributeError:
                pass

            # Perform the response correction
            ad /= ad[0].hdr['EXPTIME']
            ad /= sens_func_ad
            # Make the relevant header update
            ad.hdr['BUNIT'] = bunit

            # Now that we've made the correction, remove the superfluous
            # extra dimension from sens_func_ad and write out, if requested
            if params['write_result']:
                sens_func_ad[0].data = sens_func_ad[0].data[:, :, 0]
                try:
                    del sens_func_ad[0].WGT
                except AttributeError:
                    pass
                sens_func_ad.write(overwrite=True)
                ad.phu.set("SENSFUNC", os.path.abspath(sens_func_ad.path),
                           self.keyword_comments['SENSFUNC'])
            # sens_func_ad.reset(sens_func_regrid,
            # variance=sens_func_regrid_var)

            # Timestamp & suffix updates
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

            # import pdb;
            # pdb.set_trace()

        return adinputs

    def standardizeStructure(self, adinputs=None, **params):
        """
        The Gemini-level version of this primitive
        will try to attach an MDF because a GHOST image is
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
        Tile GHOST data into a single frame.

        This primitive will tile the SCI frames of the input images, along
        with the VAR and DQ frames if they exist.

        This primitive takes no additional parameters.
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
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by tileArrays".
                            format(ad.filename))
                adoutputs.append(ad)
                continue

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

##############################################################################
# Below are the helper functions for the user level functions in this module #
##############################################################################

    def _get_polyfit_filename(self, ad, caltype):
        """
        Gets the filename of the relevant initial polyfit file for this
        input GHOST science image

        Returns
        -------
        str/None:
            Filename (including path) of the required polyfit file
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

    def _request_bracket_arc(self, ad, before=None):
        """
        Request the 'before' or 'after' arc for the passed ad object.

        This helper function works by doing the following:
        - Append a special header keyword, 'ARCBEFOR', to the PHU. This keyword
          will be True if a before arc is requested, or False if an after arc
          is wanted.
        - getProcessedArc is the invoked, followed by the _get_cal call. The
          arc calibration association rules will see the special descriptor
          related to the 'ARCBEFOR' header keyword, and fetch an arc
          accordingly.
        - The header keyword is then deleted from the ad object, returning it
          to its original state.

        Parameters
        ----------
        before : bool
            Denotes whether to ask for the most recent arc before (True) or
            after (False) the input AD was taken. Defaults to None, at which
            point an error will be thrown.

        Returns
        -------
        arc_ad : astrodata.AstroData instance (or None)
            The requested arc. Will return None if no suitable arc is found.
        """
        if before is None:
            raise ValueError('_request_bracket_arc requires that the before '
                             'kwarg is either True or False. If you wish to '
                             'do a "standard" arc calibration fetch, simply '
                             'use getProcessedArc directly.')

        ad.phu['ARCBEFOR'] = before
        self.getProcessedArc(ad, howmany=None, refresh=True)
        arc_ad = self._get_cal(ad, 'processed_arc', )
        del ad.phu['ARCBEFOR']
        return arc_ad

    @staticmethod
    def _interp_spect(orig_data, orig_wavl, new_wavl,
                      interp='linear'):
        """
        'Re-grid' a one-dimensional input spectrum by performing simple
        interpolation on the data

        Parameters
        ----------
        orig_data : 1D numpy array or list
            The original spectrum data
        orig_wavl : 1D numpy array or list
            The corresponding wavelength values for the original spectrum data
        new_wavl : 1D numpy array or list
            The new wavelength values to re-grid the spectrum data to
        interp : str
            The interpolation method to be used. Defaults to 'linear'. Will
            accept any valid value of the ``kind`` argument to
            :any:`scipy.interpolate.interp1d`.

        Returns
        -------
        regrid_data : 1D numpy array
            The spectrum re-gridded onto the new_wavl wavelength points.
            Will have the same shape as new_wavl.
        """
        # Input checking
        orig_data = np.asarray(orig_data, dtype=orig_data.dtype)
        orig_wavl = np.asarray(orig_data, dtype=orig_wavl.dtype)
        new_wavl = np.asarray(orig_data, dtype=new_wavl.dtype)

        if orig_data.shape != orig_wavl.shape:
            raise ValueError('_interp_spect received data and wavelength '
                             'arrays of different shapes')

        interp_func = interpolate.interp1d(
            orig_wavl,
            orig_data,
            kind=interp,
            fill_value=np.nan,
            bounds_error=False,
        )
        regrid_data = interp_func(new_wavl)
        # regrid_data = np.interp(new_wavl, orig_wavl, orig_data, )
        return regrid_data

    @staticmethod
    def _regrid_spect(orig_data, orig_wavl, new_wavl,
                      waveunits='angstrom'):
        """
        Re-grid a one-dimensional input spectrum so as to conserve total flux.

        This is a more robust procedure than :meth:`_interp_spect`, and is
        designed for data with a wavelength dependence in the data units
        (e.g. erg/cm^2/s/A or similar).

        This function has been adapted from:
        http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/

        Parameters
        ----------
        orig_data : 1D numpy array or list
            The original spectrum data
        orig_wavl : 1D numpy array or list
            The corresponding wavelength values for the original spectrum data
        new_wavl : 1D numpy array or list
            The new wavelength values to re-grid the spectrum data to
        waveunits : str
            The units of the wavelength scale. Defaults to 'angstrom'.

        Returns
        -------
        egrid_data : 1D numpy array
            The spectrum re-gridded onto the new_wavl wavelength points.
            Will have the same shape as new_wavl.
        """

        spec = spectrum.ArraySourceSpectrum(wave=orig_wavl, flux=orig_data)
        f = np.ones(orig_wavl.shape)
        filt = spectrum.ArraySpectralElement(orig_wavl, f, waveunits=waveunits)
        obs = observation.Observation(spec, filt, binset=new_wavl,
                                      force='taper')
        return obs.binflux
