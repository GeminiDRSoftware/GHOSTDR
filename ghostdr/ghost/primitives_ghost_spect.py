#
#                                                                  gemini_python
#
#                                                      primitives_ghost_spect.py
# ------------------------------------------------------------------------------
import os
import numpy as np
from copy import deepcopy
import scipy
import functools
from datetime import datetime
import re

import astrodata

from geminidr.gemini.lookups import DQ_definitions as DQ

from gempy.gemini import gemini_tools as gt
from gempy.mosaic.mosaicAD import MosaicAD

from .polyfit import GhostArm, Extractor, SlitView

from .primitives_ghost import GHOST, filename_updater
from .parameters_ghost_spect import ParametersGHOSTSpect

from .lookups import polyfit_dict, line_list

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

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
                ad.write(clobber=True)

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
                dark.write(clobber=True)
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
                ad.write(clobber=True)

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

            # Convolve the flat field with the slit profile
            flat_conv = ghost_arm.slit_flat_convolve(ad[0].data,
                slit_profile=slitview.slit_profile(arm=arm),
                spatpars=spatpars[0].data, microns_pix=slitview.microns_pix,
                xpars=xpars[0].data)

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
                flatprof_ad.write(clobber=True)

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

            gt.mark_history(ad_mos, primname=self.myself(), keyword=timestamp_key)
            ad_mos.update_filename(suffix=params["suffix"], strip=True)
            adoutputs.append(ad_mos)
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

