
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
# ------------------------------------------------------------------------------
@parameter_override
class GHOSTBundle(GHOST):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the GHOSTBundle level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set(["GEMINI", "GHOST", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(GHOSTBundle, self).__init__(adinputs, **kwargs)
        self.parameters = ParametersGHOSTBundle

    def splitBundle(self, adinputs=None, **params):
        """
        This primitive breaks up a GHOST observation bundle into 3 files, one
        containing the Red camera frame, one containing the Blue camera frame,
        and another containing the Slit Viewer (SV) frames.  The Red and Blue
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

            # BEWARE: this doesn't preserve prior behaviour of putting all SV
            # exposures into a single file: test further with reduction system
            key = lambda x: (x.hdr['CAMERA'], x.hdr['EXPID'])
            extns = sorted(ad, key=key)
            for k, g in itertools.groupby(extns, key=key):
                n = astrodata.create(copy.deepcopy(ad.header[0]))
                for kw in ['NEXTEND', 'NREDEXP', 'NBLUEEXP', 'NSLITEXP']:
                    del n.phu[kw]
                for x in g: n.append(x)
                for kw in ['CAMERA', 'CCDNAME']:
                    n.phu[kw] = n[0].hdr[kw]
                n.filename = ad.filename
                binning = 'x'.join(n[0].hdr['CCDSUM'].split())
                n.update_filename(suffix='_'+binning+'_'+k[0].lower()+str(k[1]))
                log.stdinfo("   Writing {}".format(n.filename))
                n.write(clobber=True)  # should we always overwrite?

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=params['suffix'], strip=True)

        # returning [] avoids writing copy of bundle out to CWD, which then begs
        # the question: why then bother updating the timestamp & filename?
        return []#adinputs

