#
#                                                                  gemini_python
#
#                                                       primitives_ghost_slit.py
# ------------------------------------------------------------------------------
import numpy as np
from copy import deepcopy
from datetime import datetime, timedelta

from gempy.gemini import gemini_tools as gt

from .polyfit import SlitView

from .primitives_ghost import GHOST
from .primitives_ghost import filename_updater
from .parameters_ghost_slit import ParametersGHOSTSlit

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------
@parameter_override
class GHOSTSlit(GHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOSTSlit level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set(["GEMINI", "GHOST", "SLITV"])

    def __init__(self, adinputs, **kwargs):
        super(GHOSTSlit, self).__init__(adinputs, **kwargs)
        self.parameters = ParametersGHOSTSlit

    def CRCorrect(self, adinputs=None, **params):
        """
        This primitive replaces CR-affected pixels in each individual slit
        viewer image (taken from the current stream) with their equivalents
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
                # myresid.update_filename(suffix='_resid')
                # myresid.write()

                # # uncomment to output the indices for debugging
                # myindex = deepcopy(ad)
                # myindex[0].data = indices.astype(int)
                # myindex.update_filename(suffix='_index')
                # myindex.write()

                # Output and record number of pixels replaced
                nreplaced = indices.sum()
                log.stdinfo("   {}:{} : nPixReplaced = {:6d}, flux = {:.1f} -> {:.1f}"
                            "".format(ad.filename, ext.hdr['EXTVER'], nreplaced,
                                      flux_before, flux_after))
                ext.hdr['CRPIXREJ'] = (nreplaced, '# of CR pixels replaced by median')

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

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
                slitflat = gt.clip_auxiliary_data(ad, slitflat, aux_type='cal')
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
            ad.update_filename(suffix=params["suffix"], strip=True)
        return adinputs

    def stackFrames(self, adinputs=None, **params):
        """
        Combines all the extensions in a slit-viewer frame(s) into a single-
        extension AD instance.
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

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
                filename = filename_updater(ad, suffix='{:04d}'.format(index))
                adext.filename = filename
                adext.phu['ORIGNAME'] = filename
                extinputs.append(adext)

        adout = super(GHOSTSlit, self).stackFrames(extinputs, **params)[0]
        if first_filename:
            adout.phu['ORIGNAME'] = first_filename
        gt.mark_history(adout, primname=self.myself(), keyword=timestamp_key)
        adoutputs.append(adout)
        return adoutputs

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
    return reduce(lambda x, y: x + y, [np.sum(z) for z in reds + blues])

