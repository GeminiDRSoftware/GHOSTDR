
#                                                                  gemini_python
#
#                                                     primitives_ghost_bundle.py
# ------------------------------------------------------------------------------
from .primitives_ghost import GHOST, filename_updater
from .parameters_ghost_bundle import ParametersGHOSTBundle

from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

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
        self.parameters = ParametersGHOSTBundle

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
            # but one per RED/BLUE exposure which contains all SV exposures that
            # overlap with the RED/BLUE one in time (check with Jon)
            extns = [x for x in ad if x.hdr['CAMERA'].lower() == 'slit']
            _write_newfile(extns, '_slit', ad, log)

            # now do non-slitv extensions
            extns = [x for x in ad if x.hdr['CAMERA'].lower() != 'slit']
            key = lambda x: '_'+x.hdr['CAMERA'].lower()+str(x.hdr['EXPID'])
            extns = sorted(extns, key=key)
            for k, g in itertools.groupby(extns, key=key):
                _write_newfile(list(g), k, ad, log)

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=True)

        # returning [] avoids writing copy of bundle out to CWD, which then begs
        # the question: why then bother updating the timestamp & filename?
        return []#adinputs

##############################################################################
# Below are the helper functions for the primitives in this module           #
##############################################################################
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
    n = astrodata.create(copy.deepcopy(base.header[0]))
    for kw in ['NEXTEND', 'NREDEXP', 'NBLUEEXP', 'NSLITEXP']:
        del n.phu[kw]
    for x in extns: n.append(x)
    for kw in ['CAMERA', 'CCDNAME']:
        n.phu[kw] = n[0].hdr[kw]
    n.filename = base.filename
    binning = '_' + 'x'.join(n[0].hdr['CCDSUM'].split())
    n.update_filename(suffix=binning+suffix)
    log.stdinfo("   Writing {}".format(n.filename))
    n.write(overwrite=True)  # should we always overwrite?

