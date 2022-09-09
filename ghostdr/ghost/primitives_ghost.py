#
#                                                                  gemini_python
#
#                                                            primitives_ghost.py
# ------------------------------------------------------------------------------
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core.primitives_ccd import CCD
from .primitives_calibdb_ghost import CalibDBGHOST

from . import parameters_ghost

from .lookups import timestamp_keywords as ghost_stamps

from recipe_system.utils.decorators import parameter_override

import re
import astrodata
import numpy as np
# ------------------------------------------------------------------------------
_HDR_SIZE_REGEX = re.compile(r'^\[(?P<x1>[0-9]*)\:'
                             r'(?P<x2>[0-9]*),'
                             r'(?P<y1>[0-9]*)\:'
                             r'(?P<y2>[0-9]*)\]$')


# For development, let's patch in our calibration_ghost implementation to make life easier.
# Once it is finalized, we can migrate it back into the Gemini CalMgr code directly, or
# implement some more formal override option
from gemini_calmgr.cal import inst_class
import ghostdr
inst_class["GHOST"] = ghostdr.ghost.calibration_ghost.CalibrationGHOST


def filename_updater(ad, **kwargs):
    origname = ad.filename
    ad.update_filename(**kwargs)
    rv = ad.filename
    ad.filename = origname
    return rv


@parameter_override
class GHOST(Gemini, CCD, CalibDBGHOST):
    """
    Top-level primitives for handling GHOST data

    The primitives in this class are applicable to all flavours of GHOST data.
    All other GHOST primitive classes inherit from this class.
    """
    tagset = set()  # Cannot be assigned as a class

    def __init__(self, adinputs, **kwargs):
        super(GHOST, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'ghostdr.ghost.lookups'
        self._param_update(parameters_ghost)
        # Add GHOST-specific timestamp keywords
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)

    @staticmethod
    def _has_valid_extensions(ad):
        return len(ad) > 0

    def _rebin_ghost_ad(self, ad, xb, yb):
        """
        Internal helper function to re-bin GHOST data.

        .. note::
            This function is *not* a primitive. It is designed to be called
            internally by public primitives.

        This function should be used for all re-binning procedures on
        AstroData objects that will be saved during reduction. It is designed
        to handle the correct adjustment of the relevant header keywords.

        Binning is done by simply adding together the data pixels in the
        input AstroData object.

        If the input ad contains a variance plane, the re-binned variance
        plane is computed by summing over the binned variance pixels in
        quadrature. If a mask plane is present, the re-binned mask plane is
        computed by sequentially bitwise_or combining the input mask pixels.

        .. note::
            This function has been included within the GHOST primitive class
            mostly so logging is consistent. Otherwise, it could be defined as
            a @staticmethod (or just exist outside the class completely).

        Parameters
        ----------
        ad : :class:`astrodata.AstroData`
            AstroData object to be re-binned. Each extension of the object
            will be rebinned separately. A :any:`ValueError` will be thrown
            if the object's extensions are found to have different binning
            modes to one another.
        xb : :obj:`int`
            x-binning
        yb : :obj:`int`
            y-binning

        Returns
        -------
        ad : :any:`astrodata.AstroData` object
            Input AstroData object, re-binned to the requested format
        """
        log = self.log

        # Input checking
        xb = int(xb)
        yb = int(yb)
        if not isinstance(ad, astrodata.AstroData):
            raise ValueError('ad is not a valid AstroData instance')
        for ext in ad:
            if ext.hdr.get('CCDSUM') != '1 1':
                raise ValueError(
                    'Cannot re-bin data that has already been binned')

        # Re-binning
        log.stdinfo('Re-binning %s' % ad.filename)
        rows = yb
        cols = xb
        for ext in ad:
            # Do the re-binning
            # Re-bin data
            binned_array = ext.data.reshape(
                int(ext.data.shape[0] / rows), rows,
                int(ext.data.shape[1] / cols), cols
            ).sum(axis=1).sum(axis=2)

            reset_kwargs = {'check': True, }

            # Re-bin variance
            # These are summed in quadrature
            if ext.variance is not None:
                binned_var = ext.variance ** 2
                binned_var = binned_var.reshape(
                    int(ext.variance.shape[0] / rows), rows,
                    int(ext.variance.shape[1] / cols), cols
                ).sum(axis=1).sum(axis=2)
                binned_var = np.sqrt(binned_var)
                reset_kwargs['variance'] = binned_var

            # Re-bin mask
            # This can't be done in an easy one-liner - numpy bitwise_and
            # is designed to combine two distinct arrays, not combine a
            # unitary array a la, e.g. sum
            if ext.mask is not None:
                reshaped_mask = ext.mask.reshape(
                    int(ext.mask.shape[0] / rows), rows,
                    int(ext.mask.shape[1] / cols), cols
                )
                binned_mask_r = reshaped_mask[:, 0, :, :]
                for i in range(1, rows):
                    binned_mask_r = np.bitwise_or(binned_mask_r,
                                                  reshaped_mask[:, i, :, :])
                binned_mask = binned_mask_r[:, :, 0]
                for j in range(1, cols):
                    binned_mask = np.bitwise_or(binned_mask,
                                                binned_mask_r[:, :, j])
                reset_kwargs['mask'] = binned_mask

            # Alter the data values (do this all together in case one of the
            # calculations above bombs
            # ext.data = binned_array
            # ext.variance = binned_var
            # ext.mask = binned_mask
            ext.reset(binned_array, **reset_kwargs)

            # Update header values
            ext.hdr.set('CCDSUM',
                        value='%d %d' % (cols, rows,),
                        comment='Re-binned to %dx%d' % (cols, rows,))

            old_datasec = ext.hdr.get('DATASEC')
            if old_datasec:
                datasec_values = _HDR_SIZE_REGEX.match(old_datasec)
                ext.hdr.set('DATASEC',
                            value='[%d:%d,%d:%d]' %
                                  (max(int(datasec_values.group('x1')) / cols,
                                       1),
                                   max(int(datasec_values.group('x2')) / cols,
                                       1),
                                   max(int(datasec_values.group('y1')) / rows,
                                       1),
                                   max(int(datasec_values.group('y2')) / rows,
                                       1),
                                   ),
                            comment='Re-binned to %dx%d' % (cols, rows,),
                            )
            old_trimsec = ext.hdr.get('TRIMSEC')
            if old_trimsec:
                trimsec_values = _HDR_SIZE_REGEX.match(old_trimsec)
                ext.hdr.set('TRIMSEC',
                            value='[%d:%d,%d:%d]' %
                                  (max(int(trimsec_values.group('x1')) / cols,
                                       1),
                                   max(int(trimsec_values.group('x2')) / cols,
                                       1),
                                   max(int(trimsec_values.group('y1')) / rows,
                                       1),
                                   max(int(trimsec_values.group('y2')) / rows,
                                       1),
                                   ),
                            comment='Re-binned to %dx%d' % (cols, rows,),
                            )

            old_ampsize = ext.hdr.get('AMPSIZE')
            if old_ampsize:
                ampsize_values = _HDR_SIZE_REGEX.match(old_ampsize)
                ext.hdr.set('AMPSIZE',
                            value='[%d:%d,%d:%d]' %
                                  (max(int(ampsize_values.group('x1')) / cols,
                                       1),
                                   max(int(ampsize_values.group('x2')) / cols,
                                       1),
                                   max(int(ampsize_values.group('y1')) / rows,
                                       1),
                                   max(int(ampsize_values.group('y2')) / rows,
                                       1),
                                   ),
                            comment='Re-binned to %dx%d' % (cols, rows,),
                            )

        log.stdinfo('Re-binning complete')

        return ad


    def validateData(self, adinputs=None, suffix=None):
        """
        This is the data validation primitive. It checks that the instrument
        matches the primitivesClass and that there are the correct number
        of extensions.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        from gempy.gemini import gemini_tools as gt

        log = self.log
        timestamp_key = self.timestamp_keys[self.myself()]
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        missing_wcs_list = []

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by validateData".
                            format(ad.filename))
                continue

            # Check that the input is appropriate for this primitivesClass
            # Only the instrument is checked
            inst_name = ad.instrument(generic=True)
            if not inst_name in self.tagset:
                prim_class_name = self.__class__.__name__
                raise OSError("Input file {} is {} data and not suitable for "
                    "{} class".format(ad.filename, inst_name, prim_class_name))

            # Report if this is an image without square binned pixels
            if 'IMAGE' in ad.tags:
                xbin = ad.detector_x_bin()
                ybin = ad.detector_y_bin()
                if xbin != ybin:
                    log.warning("Image {} is {} x {} binned data".
                                format(ad.filename, xbin, ybin))

            if self._has_valid_extensions(ad):
                log.fullinfo("The input file has been validated: {} contains "
                             "{} extension(s)".format(ad.filename, len(ad)))
            else:
                raise OSError("The {} extension(s) in {} does not match the "
                              "number of extensions expected in raw {} "
                              "data.".format(len(ad), ad.filename, inst_name))

            # Check for WCS
            missing_wcs_list.extend([f"{ad.filename} EXTVER {ext.hdr['EXTVER']}"
                                     for ext in ad if ext.wcs is None])

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=suffix, strip=True)

        # if missing_wcs_list:
        #     msg = "The following extensions did not produce a valid WCS:\n    "
        #     msg += '\n    '.join(extstr for extstr in missing_wcs_list)
        #     raise ValueError(msg+"\n")

        return adinputs
