#
#                                                                  gemini_python
#
#                                                       primitives_ghost_slit.py
# ------------------------------------------------------------------------------
import numpy as np
from datetime import datetime, timedelta
import dateutil

from gempy.gemini import gemini_tools as gt
from geminidr.gemini.lookups import DQ_definitions as DQ

from .polyfit import SlitView

from .primitives_ghost import GHOST
from .primitives_ghost import filename_updater
from . import parameters_ghost_slit
from .lookups import polyfit_lookup

from recipe_system.utils.decorators import parameter_override
from functools import reduce

import astrodata

def parse_timestr(timestr):
    """
    Parse a time string in the format %H:%M:%S with an optional trailing .%f
    """
    if '.' not in timestr:
        timestr = timestr + '.0'
    return datetime.strptime(timestr, "%H:%M:%S.%f")

# ------------------------------------------------------------------------------
@parameter_override
class GHOSTSlit(GHOST):
    """
    Primitive class for processing GHOST slit-viewer images.
    """
    tagset = set(["GEMINI", "GHOST", "SLITV"])

    def __init__(self, adinputs, **kwargs):
        super(GHOSTSlit, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_ghost_slit)

    def darkCorrect(self, adinputs=None, **params):
        """
        Dark-correct GHOST slit observations.

        This primitive only exists to allow skipping the underlying primitive
        from DRAGONS.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        skip : bool
            Set to ``True`` to skip this primitive. Defaults to ``False``.
        dark: str/list
            name(s) of the dark file(s) to be subtracted
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # Call the underlying primitive that does the real work.
        return super().darkCorrect(adinputs, params)

    def CRCorrect(self, adinputs=None, **params):
        """
        Cosmic-ray correct slit viewer images.

        This primitive replaces CR-affected pixels in each individual slit
        viewer image (taken from the current stream) with their equivalents
        from the mean of the unaffected images.

        Cosmic rays are detected via the following algorithm:

        - Images are scaled by the total image flux
        - The median and 'median absolute deviation' (:func:`_mad <_mad>`) is
          computed for each pixel across all slit viewer frames in the stream;
        - For each slit viewer frame in the stream, a pixel is replaced by the
          corresponding mean value of the unaffected images if the pixel's
          deviation from the corresponding median is greater than some
          threshold  times the median absolute deviation for that pixel.

        Total image fluxes (computed by
        :func:`_total_obj_func <_total_obj_func>`)
        before and after pixel replacement are recorded
        in the log file, but not the file header.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        ndev: float
            number of median absolute deviations (MADs) for CR identification
        max_iters: int
            maximum number of iterations for CR removal
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        ndev = params["ndev"]
        max_iters = params["max_iters"]

        for ad in adinputs:
            filename = ad.filename
            if ad.phu.get(timestamp_key):
                log.warning(f"No changes will be made to {filename}, since "
                            "it has already been processed by CRCorrect")
                continue

            res = ad.res_mode()
            ut_date = ad.ut_date()
            binning = ad.detector_x_bin()  # x and y are equal

            fluxes = np.array([_total_obj_flux(log, res, ut_date, filename,
                                               ext.data, None, binning=binning)
                               for ext in ad])

            for iter in range(max_iters):
                total_replaced = 0
                med_flux = np.median(fluxes)
                scale_factors = med_flux / fluxes

                # Make the median slit frame by scaling individual images
                ext_stack = np.ma.masked_array([ext.data * scale for ext, scale
                                                in zip(ad, scale_factors)])
                sv_med = np.ma.median(ext_stack, axis=0)
                threshold = sv_med + ndev * _mad(ext_stack, axis=0)
                ext_stack.mask = ext_stack > threshold
                replacement_values = ext_stack.mean(axis=0)
                for i, (ext, mask) in enumerate(zip(ad, ext_stack.mask)):
                    nreplaced = mask.sum()
                    ext.data[mask] = replacement_values[mask]
                    if ext.mask is None:
                        ext.mask = np.zeros_like(ext.data, dtype=DQ.datatype)
                    ext.mask[mask] |= DQ.cosmic_ray
                    new_flux = _total_obj_flux(log, res, ut_date, filename,
                                               ext.data, None, binning=binning)
                    if nreplaced > 0:
                        log.stdinfo(f"   {filename}:{ext.hdr['EXTVER']}: "
                                    f"nPixReplaced = {nreplaced:d}, flux = "
                                    f"{fluxes[i]:.1f} -> {new_flux:.1f}")
                        fluxes[i] = new_flux
                        total_replaced += nreplaced

                if total_replaced == 0:
                    break

            # Record total number of replacements
            for ext in ad:
                ext.hdr['CRPIXREJ'] = ((ext.mask & DQ.cosmic_ray).astype(bool).sum(),
                                       '# of CR pixels replaced by median')

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)

        return adinputs

    def processSlits(self, adinputs=None, **params):
        """
        Compute the appropriate weights for each slit viewer image created
        from a bundle, so as to provide the best slit profile for each
        science exposure.

        The 'slit viewer image' for each observation is a sequence of short
        exposures of the slit viewer camera, over the entire period of the
        observation. However, only some of these will overlap each science
        exposure and it is necessary to combine different subsets in order
        to produce the best processed_slit calibration for each science frame.

        ``processSlits`` determines the appropriate weights to assign to
        each individual slit exposure, for each science exposure, and adds
        these to the SCIEXP Table for ``stackFrames`` to interpret. The
        method is to take the closest slit viewer exposure at every time
        interval during the science exposure; this is identical mathematically
        to assuming the exposures are instantaneous at their midpoints at
        linearly interpolating between them, and so copes with irregular
        temporal spacing.

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

        flat_list = params.get("flat")
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

            try:
                sciexp = ad.SCIEXP
            except AttributeError:
                log.warning(f"{ad.filename} has no SCIEXP table so all slit "
                            "viewer exposures will be combined.")
                continue

            # check that the binning is equal in x and y
            if ad.detector_x_bin() != ad.detector_y_bin():
                raise ValueError("slit viewer images must have equal binning "
                                 "in x and y directions")

            # We should process even without a slitflat; AVGEPOCH will be wrong
            # if the flux varies but we'll know which slit images to combine
            if slitflat is None:
                msg = f"Unable to find slitflat calibration for {ad.filename}"
                if self.mode == "sq":
                    raise RuntimeError(msg)
                log.warning(f"{msg}; calculations may be in error")
                sv_flat = None
            else:
                sv_flat = slitflat[0].data
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

            nslitv = len(ad)
            res = ad.res_mode()
            ut_date = ad.ut_date()
            binning = ad.detector_x_bin() # x and y are equal
            sv_duration = ad.exposure_time()
            sv_starts = [dateutil.parser.parse(f"{d}T{t}") for d, t in
                         zip(ad.hdr['DATE-OBS'], ad.hdr['UTSTART'])]
            # Relevance start and end times for each slitv image
            sv_relevant = ([datetime(2001, 1, 1)] +
                           [sv_starts[i] + 0.5 * (sv_starts[i+1] - sv_starts[i] +
                                                  timedelta(seconds=sv_duration))
                            for i in range(nslitv-1)] +
                           [datetime(3001, 1, 1)])

            fluxes = [_total_obj_flux(log, res, ut_date, ad.filename,
                                      ext.data, sv_flat, binning=binning)
                      for ext in ad]

            weights = np.zeros((len(sciexp), nslitv))
            avg_epochs = []
            for i, times in enumerate(sciexp['UTSTART', 'UTEND']):
                sc_start, sc_end = [dateutil.parser.parse(t) for t in times]
                accum_weighted_time = 0
                for j, (rel_start, rel_end, flux) in enumerate(
                        zip(sv_relevant[:-1], sv_relevant[1:], fluxes)):
                    # compute overlap fraction
                    latest_start = max(sc_start, rel_start)
                    earliest_end = min(sc_end, rel_end)
                    overlap = max((earliest_end - latest_start).total_seconds(), 0)

                    if overlap > 0:
                        weights[i, j] = overlap
                        effective_time = latest_start + 0.5 * (earliest_end -
                                                               latest_start)
                        offset = (effective_time - sc_start).total_seconds()
                        accum_weighted_time += weights[i, j] * flux * offset

                sum_of_weights = (weights[i] * fluxes).sum()
                avg_epochs.append(sc_start + timedelta(
                    seconds=accum_weighted_time / sum_of_weights))
                weights[i] /= weights[i].sum()

            sciexp['AVGEPOCH'] = [t.isoformat() for t in avg_epochs]
            for i, flux in enumerate(fluxes):
                sciexp[f"ext{i}"] = weights[:, i]
                # This could be useful for debugging purposes
                ad[i].hdr['SLITFLUX'] = (flux, "Measured slit viewer flux")

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params["suffix"], strip=True)
        return adinputs

    def stackFrames(self, adinputs=None, **params):
        """
        Combines the extensions in a slit-viewer frame into one or more
        single-extension AD instances. This is a straight weighted addition
        according to information stored in each row of the SCIEXP table. If
        there is no SCIEXP table, or create_multiple_stacks is False, then
        all the extensions are combined, with a "median" option.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        operation: str (mean|median)
            combining operation ("median" is not allowed for a weighted combne)
        create_multiple_stacks: bool
            create multiple output frames based on an attached SCIEXP table?
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        operation = params["operation"]
        create_multiple = params["create_multiple_stacks"]

        adoutputs = []
        for ad in adinputs:
            # CJS: Worth doing this check, I feel
            if 'SLITV' not in ad.tags:
                log.warning(f"{ad.filename} is not a slit-viewer image. Continuing.")
                adoutputs.append(ad)
                continue

            exptime = ad.exposure_time()
            next = len(ad)
            hdr = ad[0].hdr
            on_sky = 'CAL' not in ad.tags or 'STANDARD' in ad.tags
            # Want extension number to be last axis for broadcasting
            all_data = np.moveaxis(np.array([ext.data for ext in ad]), 0, 2)
            if any(ext.mask is None for ext in ad):
                all_mask = None
            else:
                all_mask = np.moveaxis(np.array([ext.mask for ext in ad]), 0, 2)
            if any(ext.variance is None for ext in ad):
                all_var = None
            else:
                all_var = np.moveaxis(np.array([ext.variance for ext in ad]), 0, 2)

            # extract_multiple is ignored if not an on-sky observation
            if create_multiple and on_sky:
                try:
                    sciexp_table = ad.SCIEXP
                except AttributeError:
                    log.warning("create_multiple_stacks is set to True but "
                                f"there is no SCIEXP table in {ad.filename}")
                else:
                    if operation == "median":
                        log.warning("operation='median' is invalid when "
                                    "creating multiple stacks")
                    log.stdinfo(f"Processing {ad.filename}:")
                    hdr = ad[0].hdr
                    for i, row in enumerate(sciexp_table, start=1):
                        adout = astrodata.create(ad.phu)
                        scale_factors = np.array(list(row.values())[-next:])
                        use_for = row['for']
                        log.stdinfo("  Stacking extensions for spectrograph "
                                    f"exposure(s) {use_for}")
                        log.debug(f"  Scale factors for {use_for}: {scale_factors}")
                        data, var, mask = _scale_and_stack(
                            all_data, all_var, all_mask, scale_factors)
                        adout.append(astrodata.NDAstroData(
                            data=data, mask=mask, meta={'header': hdr.copy()}))
                        adout[0].nddata.variance = var  # FIXME: can instantiate in 3.1
                        # Update keywords for calibration association purposes
                        for kw, value in zip(('DATE-OBS', 'UTSTART'),
                                             row['UTSTART'].split("T")):
                            adout.phu[kw] = value
                        adout.phu.set('ORIGTEXP', exptime, "Original exposure time")
                        adout.phu[ad._keyword_for('exposure_time')] = row['exptime']
                        adout.filename = filename_updater(
                            ad, suffix=f"_{use_for.replace(',', '_')}", strip=True)
                        adout.phu['ORIGNAME'] = adout.filename
                        adoutputs.append(adout)
                    continue

            # Regular stacking of all extensions (bypassed if SCIEXP is used)
            adout = astrodata.create(ad.phu)
            data, var, mask = _scale_and_stack(
                all_data, all_var, all_mask, np.full((next,), 1./next))
            if operation == "median":  # we keep mask but modify data and var
                data = np.median(all_data, axis=2)
                var *= 0.5 * np.pi  # according to Laplace (see nddops.py)
            adout.append(astrodata.NDAstroData(
                data=data, mask=mask, meta={'header': hdr.copy()}))
            adout[0].nddata.variance = var  # FIXME: can instantiate in 3.1
            adout.phu['ORIGNAME'] = ad.phu['ORIGNAME'] or ad.filename
            adoutputs.append(adout)

        for ad in adoutputs:
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            # This stuff is in the PHU so delete vestigial info from extensions
            for kw in ('DATE-OBS', 'UTSTART', 'UTEND', 'EXPUTST', 'EXPUTEND'):
                del ad[0].hdr[kw]
        return adoutputs

##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################
def _mad(data, axis=None, keepdims=False):
    """
    Median Absolute Deviation: a "Robust" version of standard deviation.

    The median absolute deviation of a sample is the median of data devations
    from the data median:

    .. math::
        \\textrm{MAD} = \\textrm{median} ( | X_i - \\textrm{median}(X) | )

    For further details, see:
    https://en.wikipedia.org/wiki/Median_absolute_deviation

    Parameters
    ----------
    data : list or numpy array
        Data to get the 'median absolute variation' for
    axis : int or None
        Axis along which to compute the MAD. Defaults to None (i.e. MAD
        is computed across all data points).
    keepdims: bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original data.

    Returns
    -------
    float
        The MAD of the data, along the requested axis.
    """
    return np.ma.median(np.absolute(data - np.ma.median(
        data, axis=axis, keepdims=True)), axis=axis, keepdims=keepdims)


def _scale_and_stack(all_data, all_var=None, all_mask=None, scale_factors=None):
    """
    Stacks multiple data planes with predetermined scalings. The variance
    and masks are propagated if they exists: in the case of the mask, the
    bitwise_or operation is performed on all data planes with non-zero scale
    factors.

    Parameters
    ----------
    all_data: float array (ny, nx, next)
        array of data from all extensions
    all_data: variance array (ny, nx, next) or None
        array of variance from all extensions
    all_mask: uint16 array (ny, nx, next) or None
        array of masks from all extensions
    scale_factors: float array (next)
        array of scale factors to apply to each extension

    Returns
    -------
    out_data: float array (ny, nx)
        output data array
    out_var: float array (ny, nx) or None
        output variance array
    out_mask: uint16 array (ny, nx) or None
        output mask
    """
    out_data = np.sum(all_data * scale_factors, axis=-1)
    if all_var is None:
        out_var = None
    else:
        out_var = np.sum(all_var * scale_factors**2, axis=-1)
    if all_mask is None:
        out_mask = None
    else:
        out_mask = np.bitwise_or.reduce(
            all_mask & np.where(scale_factors > 0, -1, 0),
            axis=-1).astype(all_mask.dtype)
    return out_data, out_var, out_mask


def _total_obj_flux(log, res, ut_date, filename, data, flat_data=None, binning=2):
    """
    Combined red/blue object flux calculation.

    Uses the :any:`polyfit.slitview.SlitView` object to
    determine (potentially sky-subtracted) total object flux. In high-resolution
    mode, the concurrent arc profile is returned as an "object" profile,
    so we discard it explicitly from this calculation.

    Sky subtraction occurs if the ``flat_data`` parameter is not :any:`None`.

    Parameters
    ----------
    res: string
        Either ``'high'`` or ``'std'``.
    data: :class:`numpy.ndarray`
        The slit viewer image data from which to extract the object profiles
    flat_data: :class:`numpy.ndarray`/None
        The bias-/dark-corrected slit view flat field image used to determine
        sky background levels (may be ``None`` if sky subtraction not
        needed).

    Returns
    -------
    flux: float
        The object flux, summed, and potentially sky-subtracted.
    """
    sky_correction = flat_data is not None
    slitv_fn = polyfit_lookup.get_polyfit_filename(log, 'slitv',
                                                   res, ut_date, filename,
                                                   'slitvmod')
    slitvpars = astrodata.open(slitv_fn)
    svobj = SlitView(data, flat_data, slitvpars.TABLE[0], mode=res,
                     microns_pix=4.54*180/50, binning=binning)  # OK to pass None for flat
    reds = svobj.object_slit_profiles(
        'red', correct_for_sky=sky_correction, append_sky=False,
        normalise_profiles=False)
    blues = svobj.object_slit_profiles(
        'blue', correct_for_sky=sky_correction, append_sky=False,
        normalise_profiles=False)

    # discard the arc profiles if high res
    if res == 'high':
        blues = blues[:1]
        reds = reds[:1]
    return reduce(lambda x, y: x + y, [np.sum(z) for z in reds + blues])

