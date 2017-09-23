#
#                                                                  gemini_python
#
#                                                            primitives_ghost.py
# ------------------------------------------------------------------------------
import os
import numpy as np
from copy import deepcopy
import scipy
import functools
from datetime import datetime, timedelta

import astrodata

from geminidr.gemini.lookups import DQ_definitions as DQ

from gempy.gemini import gemini_tools as gt
from gempy.mosaic.mosaicAD import MosaicAD

from .polyfit import GhostArm, Extractor, SlitView

from geminidr.core import CCD
from geminidr.gemini.primitives_gemini import Gemini
from .primitives_calibdb_ghost import CalibDBGHOST

from .parameters_ghost import ParametersGHOST

from .lookups import timestamp_keywords as ghost_stamps, polyfit_dict, line_list

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------
@parameter_override
class GHOST(Gemini, CCD, CalibDBGHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOST level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set(["GEMINI", "GHOST"])

    def __init__(self, adinputs, **kwargs):
        super(GHOST, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'ghostdr.ghost.lookups'
        self.parameters = ParametersGHOST
        # Add GHOST-specific timestamp keywords
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)

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
            ad.filename = gt.filename_updater(ad, suffix=params["suffix"],
                                              strip=True)
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

    def correctSlitCosmics(self, adinputs=None, **params):
        """
        This primitive replaces CR-affected pixels in each individual slit
        viewer image (taken from the current stream) with their equivalents\
        from the median frame of those images.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Check that all SLITV inputs are single-extension images
        #if any([('SLITV' in ad.tags and len(ad)!=1) for ad in adinputs]):
        #    raise IOError("All input SLITV images must have a single extension")

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by correctSlitCosmics".
                            format(ad.filename))
                continue

            if 'SLITV' not in ad.tags:
                log.warning("No changes will be made to {}, since it is not "
                            "a slit viewer frame".format(ad.filename))
                continue

            # Make the median slit frame.
            # ext_stack = np.array([ad[0].data for ad in adinputs])
            ext_stack = np.array([ext.data for ext in ad])
            sv_med = np.median(ext_stack, axis=0)
            sigma = _mad(ext_stack, axis=0) * 20
            res = ad.res_mode()

            for ext in ad:
                addata = ext.data

                # pre-CR-corrected flux computation
                flux_before = _total_obj_flux(res, addata, None)

                # replace CR-affected pixels with those from the median slit
                # viewer image (subtract the median slit frame and threshold
                # against the residuals); not sure VAR/DQ planes appropriately
                # handled here
                residuals = abs(addata - sv_med)
                indices = residuals > sigma
                addata[indices] = sv_med[indices]

                # post-CR-corrected flux computation
                flux_after = _total_obj_flux(res, addata, None)

                # CJS: since addata is a reference, the CR fix changes the data
                # in place and there's no need to reassign ad[0].data = addata

                # # uncomment to output the residuals for debugging
                # myresid = deepcopy(ad)
                # myresid[0].data = residuals
                # myresid.filename = gt.filename_updater(ad, suffix='_resid')
                # myresid.write()

                # # uncomment to output the indices for debugging
                # myindex = deepcopy(ad)
                # myindex[0].data = indices.astype(int)
                # myindex.filename = gt.filename_updater(ad, suffix='_index')
                # myindex.write()

                # Output and record number of pixels replaced
                nreplaced = indices.sum()
                log.stdinfo("   {}:{} : nPixReplaced = {:6d}, flux = {:.1f} -> {:.1f}"
                            "".format(ad.filename, ext.hdr['EXTVER'], nreplaced,
                                      flux_before, flux_after))
                ext.hdr['CRPIXREJ'] = (nreplaced, '# of CR pixels replaced by median')

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.filename = gt.filename_updater(ad, suffix=params["suffix"],
                                              strip=True)

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
        # Get processed slits, slitFlats, and Xmods
        # slits and slitFlats may be provided as parameters
        slit_list = params["slit"]
        if slit_list is None:
            # CJS: This populates the calibrations cache (dictionary) with
            # "processed_slit" filenames for each input AD
            self.getProcessedSlit(adinputs)
            # This then gets those filenames
            slit_list = [self._get_cal(ad, 'processed_slit')
                         for ad in adinputs]

        flat_list = params["slitflat"]
        if flat_list is None:
            self.getProcessedSlitFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        self.getProcessedXmod(adinputs)
        xmod_list = [self._get_cal(ad, 'processed_xmod') for ad in adinputs]

        # TODO: Have gt.make_lists handle multiple auxiliary lists?
        # CJS: Here we call gt.make_lists. This has only been designed to work
        # with one auxiliary list at present, hence the three calls. This
        # produces two lists of AD objects the same length, one of the input
        # ADs and one of the auxiliary files, from the list
        # of filenames (or single passed parameter). Importantly, if multiple
        # auxiliary frames are the same, then the file is opened only once and
        # the reference to this AD is re-used, saving speed and memory.
        _, slit_list = gt.make_lists(adinputs, slit_list, force_ad=True)
        _, flat_list = gt.make_lists(adinputs, flat_list, force_ad=True)
        _, xmod_list = gt.make_lists(adinputs, xmod_list, force_ad=True)

        for ad, slit, flat, xpars in zip(adinputs, slit_list,
                                        flat_list, xmod_list):
            log.info(ad.info())

            # CJS: failure to find a suitable auxiliary file (either because
            # there's no calibration, or it's missing) places a None in the
            # list, allowing a graceful continuation.
            if slit is None or flat is None or xpars is None:
                log.warning("Unable to find calibrations for {}; "
                            "skipping".format(ad.filename))
                continue

            # CJS: Changed to log.debug() and changed the output
            log.debug("Slit parameters: ")
            log.debug("   processed_list: {}".format(slit.filename))
            log.debug("   processed_flat: {}".format(flat.filename))

            res_mode = ad.res_mode()
            arm = GhostArm(arm=ad.arm(), mode=res_mode)

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

            arm.spectral_format_with_matrix(xpars[0].data, wpars[0].data,
                        spatpars[0].data, specpars[0].data, rotpars[0].data)
            sview = SlitView(slit[0].data, flat[0].data, mode=res_mode)

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
            ad.filename = gt.filename_updater(ad, suffix=params["suffix"],
                                              strip=True)
            if params["write_result"]:
                ad.write(clobber=True)

        return adinputs

    def findApertures(self, adinputs=None, **params):
        """
        Locate the apertures within a GHOST frame, and write out polyfit-
        compliant FITS files to the calibrations system
        
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

        adoutputs = []
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

            # Minimal PrimaryHDU, but AD must recognize as "GHOST"
            ad_xmod = astrodata.create({'INSTRUME': 'GHOST'})
            ad_xmod.append(fitted_params)
            ad_xmod.filename = ad.filename
            # CJS: Need to add data_label for storeAsCalibration()
            ad_xmod.phu['DATALAB'] = ad.data_label()
            gt.mark_history(ad_xmod, primname=self.myself(), keyword=timestamp_key)

            adoutputs.append(ad_xmod)
        return adoutputs

    def fitWavelength(self, adinputs=None, **params):
        """
        Fit wavelength solution to a GHOST ARC frame
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        self.getProcessedXmod(adinputs)
        xmod_list = [self._get_cal(ad, 'processed_xmod') for ad in adinputs]

        adoutputs = []
        for ad, xpars in zip(*gt.make_lists(adinputs, xmod_list, force_ad=True)):
            if not {'GHOST', 'PROCESSED', 'ARC'}.issubset(ad.tags):
                log.warning("fitWavelength is only run on prepared GHOST arc"
                            " files - skipping {}".format(ad.filename))
                continue

            if self.timestamp_keys["extractProfile"] not in ad.phu.keywords:
                log.warning("extractProfile has not been run on {} - "
                            "skipping".format(ad.filename))
                continue

            if xpars is None:
                log.warning("Could not find processed_xmod calibration for "
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
            arm.spectral_format_with_matrix(xpars[0].data,
                                            wpars[0].data,
                                            spatpars[0].data,
                                            specpars[0].data,
                                            rotpars[0].data)

            extractor = Extractor(arm, None)  # slitview=None for this usage
            lines_out = extractor.find_lines(ad[0].data, arcwaves, inspect=False)

            #import pdb; pdb.set_trace()
            fitted_params, wave_and_resid = arm.read_lines_and_fit(
                wpars[0].data, lines_out)

            # Much like the solution for findApertures, create a minimum-spec
            # AstroData object to prepare the result for storage in the
            # calibrations system
            # CJS: Renamed to ad_wfit from ad_xmod
            ad_wfit = astrodata.create({'INSTRUME': 'GHOST'})
            ad_wfit.append(fitted_params)
            ad_wfit.filename = ad.filename
            # CJS: Need to add data_label for storeAsCalibration()
            ad_wfit.phu['DATALAB'] = ad.data_label()
            gt.mark_history(ad_wfit, primname=self.myself(), keyword=timestamp_key)

            adoutputs.append(ad_wfit)
        return adoutputs

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

        flat_list = params["slitflat"]
        if flat_list is None:
            self.getProcessedSlitFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        objflat_list = params["flat"]
        if objflat_list is None:
            self.getProcessedFlat(adinputs)
            objflat_list = [self._get_cal(ad, 'processed_flat')
                         for ad in adinputs]

        self.getProcessedXmod(adinputs)
        xmod_list = [self._get_cal(ad, 'processed_xmod') for ad in adinputs]

        # TODO: Have gt.make_lists handle multiple auxiliary lists?
        _, slit_list = gt.make_lists(adinputs, slit_list, force_ad=True)
        _, flat_list = gt.make_lists(adinputs, flat_list, force_ad=True)
        _, objflat_list = gt.make_lists(adinputs, objflat_list, force_ad=True)
        _, xmod_list = gt.make_lists(adinputs, xmod_list, force_ad=True)

        for ad, slit, flat, obj_flat, xpars in zip(adinputs, slit_list,
                                        flat_list, objflat_list, xmod_list):
            log.info(ad.info())

            # CJS: failure to find a suitable auxiliary file (either because
            # there's no calibration, or it's missing) places a None in the
            # list, allowing a graceful continuation.
            if slit is None or flat is None or xpars is None or obj_flat is None:
                log.warning("Unable to find calibrations for {}; "
                            "skipping".format(ad.filename))
                continue

            try:
                wpars = astrodata.open(self._get_polyfit_filename(ad, 'wavemod'))
                spatpars = astrodata.open(self._get_polyfit_filename(ad, 'spatmod'))
                specpars = astrodata.open(self._get_polyfit_filename(ad, 'specmod'))
                rotpars = astrodata.open(self._get_polyfit_filename(ad, 'rotmod'))
            except IOError:
                log.warning("Cannot open required initial model files for {};"
                            " skipping".format(ad.filename))
                continue

            res_mode = ad.res_mode()
            arm = GhostArm(arm=ad.arm(), mode=res_mode)
            arm.spectral_format_with_matrix(xpars[0].data,
                                            wpars[0].data,
                                            spatpars[0].data,
                                            specpars[0].data,
                                            rotpars[0].data)

            sview = SlitView(slit[0].data, flat[0].data, mode=res_mode)

            extractor = Extractor(arm, sview)
            print obj_flat.filename, obj_flat[0].data.shape
            extracted_flux, extracted_var = extractor.two_d_extract(
                obj_flat[0].data, extraction_weights=ad[0].WGT)
            # import pdb; pdb.set_trace()

            # Normalised extracted flat profile
            med = np.median(extracted_flux)
            extracted_flux /= med
            extracted_var /= med**2

            flatprof_ad = deepcopy(ad)
            flatprof_ad.filename = gt.filename_updater(flatprof_ad,
                            suffix='_extractedFlatProfile, strip=True')
            flatprof_ad[0].reset(extracted_flux, mask=None,
                                 variance=extracted_var)
            if params["write_result"]:
                flatprof_ad.write(clobber=True)

            # Divide the flat field through the science data
            # Arithmetic propagates VAR correctly
            ad /= flatprof_ad

        return adinputs

    def processSlits(self, adinputs=None, **params):
        """
        This primitive computes the mean exposure epoch for an input SLITV
        image (time series of slit-viewer images) and writes it into the PHU

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        slitflat: str/None
            name of the slitflat to use (if None, use the calibration
            system)
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        flat_list = params["flat"]
        if flat_list is None:
            self.getProcessedSlitFlat(adinputs)
            flat_list = [self._get_cal(ad, 'processed_slitflat')
                         for ad in adinputs]

        for ad, slitflat in zip(*gt.make_lists(adinputs, flat_list,
                                               force_ad=True)):
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by processSlits".
                            format(ad.filename))
                continue

            if slitflat is None:
                log.warning("Unable to find slitflat calibration for {}; "
                            "skipping".format(ad.filename))
                continue
            else:
                sv_flat = slitflat[0].data

            # accumulators for computing the mean epoch
            sum_of_weights = 0.0
            accum_weighted_time = 0.0

            # Check the inputs have matching binning and SCI shapes.
            try:
                gt.check_inputs_match(adinput1=ad, adinput2=slitflat,
                                      check_filter=False)
            except ValueError:
                # This is most likely because the science frame has multiple
                # extensions and the slitflat needs to be copied
                slitflat = gt.clip_auxiliary_data(ad, slitflat, aux_type='cal',
                                    keyword_comments=self.keyword_comments)
                # An Error will be raised if they don't match now
                gt.check_inputs_match(ad, slitflat, check_filter=False)

            # get science start/end times
            sc_start = datetime.strptime(ad.phu['UTSTART'], "%H:%M:%S.%f")
            sc_end = datetime.strptime(ad.phu['UTEND'], "%H:%M:%S.%f")

            res = ad.res_mode()
            for ext in ad:
                sv_start = datetime.strptime(ext.hdr['EXPUTST'], "%H:%M:%S.%f")
                sv_end = datetime.strptime(ext.hdr['EXPUTEND'], "%H:%M:%S.%f")

                # compute overlap percentage and slit view image duration
                latest_start = max(sc_start, sv_start)
                earliest_end = min(sc_end, sv_end)
                overlap = (earliest_end - latest_start).seconds
                overlap = 0.0 if overlap < 0.0 else overlap  # no overlap edge case
                sv_duration = ext.hdr['EXPTIME']
                overlap /= sv_duration  # convert into a percentage

                # compute the offset (the value to be weighted), in seconds,
                # from the start of the science exposure
                offset = 42.0  # init value: overridden if overlap, else 0-scaled
                if sc_start <= sv_start and sv_end <= sc_end:
                    offset = (sv_start - sc_start).seconds + sv_duration / 2.0
                elif sv_start < sc_start:
                    offset = overlap * sv_duration / 2.0
                elif sv_end > sc_end:
                    offset = overlap * sv_duration / 2.0
                    offset += (sv_start - sc_start).seconds

                # add flux-weighted offset (plus weight itself) to accumulators
                flux = _total_obj_flux(res, ext.data, sv_flat)
                weight = flux * overlap
                sum_of_weights += weight
                accum_weighted_time += weight * offset

            # final mean exposure epoch computation
            if sum_of_weights > 0.0:
                mean_offset = accum_weighted_time / sum_of_weights
                mean_offset = timedelta(seconds=mean_offset)
                # write the mean exposure epoch into the PHU
                sc_start = datetime.strptime(ad.phu['UTSTART'], "%H:%M:%S.%f")
                mean_epoch = sc_start + mean_offset
                ad.phu['AVGEPOCH'] = (  # hope this keyword string is ok
                    mean_epoch.strftime("%H:%M:%S.%f")[:-3],
                    'Mean Exposure Epoch')

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.filename = gt.filename_updater(ad, suffix=params["suffix"],
                                              strip=True)
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
            ad.filename = gt.filename_updater(ad, suffix=params["suffix"],
                                              strip=True)
        return adinputs

    def stackSlitFrames(self, adinputs=None, **params):
        """
        Combines all the extensions in a slit-viewer frame(s) into a single-
        extension AD instance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        operation = params["operation"]
        reject_method = params["reject_method"]

        # Keep hold of the first SLITV input's filename so the stacked
        # output has a sensible name.
        first_filename = None

        # CJS: This could be rewritten to produce one output AD for each
        # input AD by combining the extensions from each input separately.
        # That behaviour would not need the addToList/getList primitives.
        # But this is a better match to how the original code behaved.
        adoutputs = []
        extinputs = []
        for ad in adinputs:
            # CJS: Worth doing this check, I feel
            if 'SLITV' not in ad.tags:
                log.warning("{} is not a slit-viewer image. Continuing.".
                            format(ad.filename))
                adoutputs.append(ad)
                continue

            if not first_filename:
                first_filename = ad.phu['ORIGNAME'] or ad.filename
            # DQ plane is still needed so call stackFrames for ease
            # CJS: This is ugly but should go with pythonic stacking
            for index, ext in enumerate(ad, start=1):
                adext = deepcopy(ext)
                filename = gt.filename_updater(ad, suffix='{:04d}'.format(index))
                adext.filename = filename
                adext.phu['ORIGNAME'] = filename
                extinputs.append(adext)

        # CJS: Could simply pass parameters as **params here
        adout = self.stackFrames(extinputs, operation=operation,
                                 reject_method=reject_method)[0]
        if first_filename:
            adout.phu['ORIGNAME'] = first_filename
        gt.mark_history(adout, primname=self.myself(), keyword=timestamp_key)
        adoutputs.append(adout)
        return adoutputs

    def standardizeStructure(self, adinputs=None, **params):
        """
        CJS: Only exists now to add DATASEC keyword to SLITV frames, which
        is missing in the simulated data.
        No longer promotes extensions of SLITV images to full AD instances
        since stackSlitFrames() handles this.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        for ad in adinputs:
            if 'SLITV' in ad.tags:
                ad.hdr['DATASEC'] = '[1:{1},1:{0}]'.format(*ad[0].data.shape)
        return adinputs

    # CJS: Primitive has been renamed for consistency with other instruments
    def tileArrays(self, adinputs=None, **params):
        """
        This primitive will tile the SCI frames of the input images, along
        with the VAR and DQ frames if they exist.
        """
        def simple_mosaic_function(ad):
            """
            This should probably go into MosaicAD as the default function.
            Being discussed within the team.
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
            ad_mos.filename = gt.filename_updater(ad_mos, suffix=params["suffix"],
                                              strip=True)
            adoutputs.append(ad_mos)
        return adoutputs

    # validateData() removed since inherited Standardize method will handle it

    def _get_polyfit_filename(self, ad=None, cal_type=None):
        """
        Gets the filename of the relevant initial polyfit file for this
        input GHOST science image

        Returns
        -------
        str/None: Filename of the required polyfit file
        """
        log = self.log
        polyfit_dir = os.path.join(os.path.dirname(polyfit_dict.__file__),
                                   'Polyfit')

        # CJS: This is a method that only exists *if* the input is of type
        # GHOST, so no need to check
        xbin = ad.detector_x_bin()
        ybin = ad.detector_y_bin()
        arm = ad.arm()
        res_mode = ad.res_mode()
        key = 'GHOST_{}_{}_{}_{}'.format(xbin, ybin, arm, res_mode)

        try:
            poly_dict = getattr(polyfit_dict, '{}_dict'.format(cal_type))
        except AttributeError:
            log.warning("Invalid polyfit calibration type ({}) requested for "
                        "{}".format(cal_type, ad.filename))
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

##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################
def _mad(data, axis=None):
    """
    Median Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    return np.median(np.absolute(data - np.median(data, axis)), axis)

def _total_obj_flux(res, data, flat_data=None):
    """
    combined red/blue object flux calculation. uses the slitview object to
    determine (potentially sky-subtracted) object profiles. in high res
    mode, the arc profile is returned as an "object" profile, so we discard
    it explicitly from this calculation
    
    Parameters
    ----------
    res: string
        either 'high' or 'std'
    data: np.ndarray
        the slit viewer image data from which to extract the object profiles
    flat_data: np.ndarray/None
        the bias-/dark-corrected slit view flat field image used to de-
        termine sky background levels (may be None if sky subtraction not
        needed)

    Returns
    -------
    flux: float
        the object flux, summed, and potentially sky-subtracted
    """
    sky_correction = flat_data is not None
    svobj = SlitView(data, flat_data, mode=res)  # OK to pass None for flat
    reds = svobj.object_slit_profiles('red',  # noqa
        correct_for_sky=sky_correction, append_sky=False,
                                      normalise_profiles=False)
    blues = svobj.object_slit_profiles('blue',  # noqa
        correct_for_sky=sky_correction, append_sky=False,
                                      normalise_profiles=False)

    # discard the arc profiles if high res
    if res == 'high':
        blues = blues[:1]
        reds = reds[:1]
    return reduce(lambda x, y: x+y, [np.sum(z) for z in reds+blues])
