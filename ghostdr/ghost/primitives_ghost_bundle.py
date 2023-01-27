
#                                                                  gemini_python
#
#                                                     primitives_ghost_bundle.py
# ------------------------------------------------------------------------------
from .primitives_ghost import GHOST, filename_updater
from . import parameters_ghost_bundle

from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from collections import Counter
import astrodata
import copy
import itertools
from astropy.io.fits import PrimaryHDU, Header

# ------------------------------------------------------------------------------
@parameter_override
class GHOSTBundle(GHOST):
    """
    Primitives for unpacking GHOST observation bundle files.
    """
    tagset = set(["GEMINI", "GHOST", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(GHOSTBundle, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_ghost_bundle)

    def splitBundle(self, adinputs=None, **params):
        """
        Break a GHOST observation bundle into individual exposures.

        This primitive breaks up a GHOST observation bundle into 3 files: one
        containing the Red camera frame, one containing the Blue camera frame,
        and another containing the Slit Viewer (SV) frames.

        The Red and Blue
        output files are MEF because each amp quadrant is in its own extension,
        while the SV output file will contain all SV exposures taken during the
        observation run and will thus be single-extension for zero-duration
        BIAS observations, but it may also be a MEF for other observation types
        due to their longer exposures.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by {}".format(
                            ad.filename, self.myself()))
                continue

            log.stdinfo("Unbundling {}:".format(ad.filename))

            # as a special case, write all slitv extns to a single file
            # TODO: may need to make multiple SV files, not one per SV exposure
            # FIXME: better way to detect slit exposures than by camera?
            # but one per RED/BLUE exposure which contains all SV exposures that
            # overlap with the RED/BLUE one in time (check with Jon)
            extns = [x for x in ad if (x.arm() == 'slitv' and x.shape)]
            if len(extns) > 0:
                _write_newfile(extns, '_slit', ad, log)

            # now do non-slitv extensions
            extns = [x for x in ad if x.arm() != 'slitv']
            key = lambda x: f"_{x.hdr['CAMERA'].lower()}{x.hdr['EXPID']:03d}"
            extns = sorted(extns, key=key)
            for k, g in itertools.groupby(extns, key=key):
                _write_newfile(list(g), k, ad, log)

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=True)

        # returning [] avoids writing copy of bundle out to CWD, which then begs
        # the question: why then bother updating the timestamp & filename?
        return []#adinputs

    def validateData(self, adinputs=None, suffix=None):
        """
        GHOSTBundle-specific version of validateData to ignore the invalid WCS
        exception.
        """
        try:
            super().validateData(adinputs, suffix=suffix)
        except ValueError as e:
            if 'valid WCS' not in str(e):
                raise
        return adinputs


##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################
def _get_common_hdr_value(base, extns, key):
    """
    Helper function to get a common header value from a list of extensions.

    Parameters
    ----------
    base : :any:`astrodata.AstroData`
        Original AstroData instance that may contain useful headers
    extns : iterable of :any:`astrodata.Astrodata`
        AstroData extensions to be examined
    key : str
        The FITS header key to find in each extension 

    Raises
    ------
    KeyError
        If the ``key`` parameter has multiple values across the extensions
    """

    # Get the keyword from every extension
    vals = [x.hdr.get(key) for x in extns]
    c = Counter(vals)
    # Not all extensions may not contain the keyword,
    # but we don't care about blanks
    del c[None]
    # If the keyword doesn't exist at all in the extensions,
    # then use the base value instead
    if len(c) == 0:
        return base.phu.get(key)
    # Check that every extension has the same value for the keyword
    most_common = c.most_common()
    base = most_common[0]
    if len(most_common) != 1:
        # Ignore single errors in this header, as long as the most common
        # value is more than half
        if base[1] < len(vals) // 2:
          raise KeyError('multiple values for ' + key + " " + str(vals))

        for val in most_common[1:]:
          if val[1] == 1:
            print('Ignoring single error in', key, 'header')
          else:
            raise KeyError('multiple values for ' + key + " " + str(vals))

    return base[0]

def _get_hdr_values(extns, key):
    """
    Helper function to get the all header values from a list of extensions.
    The return value is a dict keyed on the EXPID value.

    Parameters
    ----------
    extns : iterable of :any:`astrodata.Astrodata`
        AstroData extensions to be examined
    key : str
        The FITS header key to find in the list of extensions
    """

    return {x.hdr.get('EXPID'): x.hdr.get(key) for x in extns}


def _write_newfile(extns, suffix, base, log):
    """
    Helper function to write sub-files out from a MEF bundle.

    Parameters
    ----------
    extns : iterable of :any:`astrodata.Astrodata`
        AstroData extensions to be appended to the new file
    suffix : str
        Suffix to be appended to file name
    base : :any:`astrodata.AstroData`
        Original AstroData instance to base the new file on. The file's
        primary header unit will be used, as will the base of the filename.
    log : AstroData logging object
        Log for recording actions. This should be the log in use in the calling
        primitive.

    Raises
    ------
    AssertionError
        If the ``extns`` parameter is :any:`None`, or empty
    """
    assert extns and len(extns) > 0
    # Start with a copy of the base PHU
    n = astrodata.create(copy.deepcopy(base.phu))

    # But also look for the extension with an empty data array,
    # because this is the real PHU from the original file before
    # it was mashed into a giant MEF.
    for x in extns:
        if (x.hdr.get('NAXIS') == 0) or (x.data.size == 0):
            phu = PrimaryHDU(data=None, header=copy.deepcopy(x.hdr))
            n = astrodata.create(phu)

    # Copy some important keywords into each separate file if they
    # aren't already there
    for kw in ['INSTRUME', 'TELESCOP', 'DATALAB', 'GEMPRGID', 'OBSID', 'UTC-OBS', 'RA', 'DEC']:
        try:
            if not n.phu.get(kw):
                n.phu[kw] = base.phu.get(kw)
        except KeyError:
            pass

    # Remove some keywords that are only relevant to the bundle
    for kw in ['NEXTEND', 'NREDEXP', 'NBLUEEXP', 'NSLITEXP']:
        try:
            del n.phu[kw]
        except KeyError:
            pass

    # Append the extensions that contain pixel data
    for x in extns:
        if (x.hdr.get('NAXIS') > 0) and (x.data.size > 0):
            n.append(x)

    # Collate headers into the new PHU
    for kw in ['CAMERA', 'CCDNAME',
               'CCDSUM', 'DETECTOR',
               'OBSTYPE', 'SMPNAME']:
        n.phu.set(kw, _get_common_hdr_value(base, extns, kw))
    vals = _get_hdr_values(extns, 'DATE-OBS')
    n.phu.set('DATE-OBS', vals[min(vals.keys())])
    vals = _get_hdr_values(extns, 'UTSTART')
    n.phu.set('UTSTART', vals[min(vals.keys())])
    vals = _get_hdr_values(extns, 'UTEND')
    n.phu.set('UTEND', vals[max(vals.keys())])

    # Construct a filename
    n.filename = base.filename
    ccdsum = n.phu.get('CCDSUM')
    binning = '_' + 'x'.join(ccdsum.split())
    n.update_filename(suffix=binning+suffix)

    # MCW 190813 - Update the ORIGNAME of the file
    # Otherwise, every time we do a 'strip' file rename, the base file name
    # will go back to being the MEF bundle file name, and things will
    # quickly start to overlap each other
    n.phu['ORIGNAME'] = n.filename

    # CJS 20221128: to ensure that processed cals from the different arms
    # have different data labels before going in the archive
    n.phu['DATALAB'] += f"-{n.phu['CAMERA']}"
    if n.phu['CAMERA'] != "SLITV":
        n.phu['DATALAB'] += f"-{suffix[-3:]}"

    log.stdinfo("   Writing {}".format(n.filename))
    n.write(overwrite=True)  # should we always overwrite?

