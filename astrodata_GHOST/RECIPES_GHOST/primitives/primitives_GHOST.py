from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors
from gempy.gemini import gemini_tools as gt
from astrodata_GHOST.ADCONFIG_GHOST.lookups import timestamp_keywords as ghost_stamps
from gempy.gemini.eti.gireduceparam import subtract_overscan_hardcoded_params
from pyraf import iraf

from primitives_GMOS import GMOSPrimitives

class GHOSTPrimitives(GMOSPrimitives):
    """
    Class containing all the GHOST primitives.  It inherits all the primitives
    from the GEMINIPrimitives class (which itself inherits a series of 
    primitives from RECIPE_Gemini/primitives.)
    """

    astrotype = "GHOST"

    def init(self, rc):
        GMOSPrimitives.init(self, rc)
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)
        return rc

    def subtractOverscan(self, rc):
        iraf.setVerbose(value=2)
        subtract_overscan_hardcoded_params['order'] = 1
        return GMOSPrimitives.subtractOverscan(self, rc)

    def standardizeHeaders(self, rc):
        """
        This primitive is used to standardize the headers of GHOST data,
        specifically.
        """

        log = logutils.get_logger(__name__)
        log.debug(gt.log_message("primitive", "standardizeHeaders", "starting"))
        #timestamp_key = self.timestamp_keys["standardizeHeaders"]
        rc.run('standardizeGeminiHeaders')
        yield rc

