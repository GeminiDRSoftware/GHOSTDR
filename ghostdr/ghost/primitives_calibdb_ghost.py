#
#                                                                  gemini_python
#
#                                                    primitives_calibdb_ghost.py
# ------------------------------------------------------------------------------
from gempy.gemini import gemini_tools as gt

from geminidr.core import CalibDB
from geminidr.core.primitives_calibdb import REQUIRED_TAG_DICT
from . import parameters_calibdb_ghost


from recipe_system.utils.decorators import parameter_override

# Need to add our special calibrator types to the REQUIRED_TAG_DICT

REQUIRED_TAG_DICT['processed_slitflat'] = ['PROCESSED', 'FLAT', 'SLITV']
REQUIRED_TAG_DICT['processed_slit'] = ['PROCESSED', 'SLITV']

# ------------------------------------------------------------------------------
@parameter_override
class CalibDBGHOST(CalibDB):
    """
    This is the class containing all of the calibration bookkeeping primitives
    for the CalibDBGHOST level of the type hierarchy tree. It inherits all
    the primitives from the level above
    """
    tagset = set()  # Not allowed to be a selected as a primitivesClass

    def __init__(self, adinputs, **kwargs):
        super(CalibDBGHOST, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'ghostdr.ghost.lookups'
        self._param_update(parameters_calibdb_ghost)

    def getProcessedArc(self, adinputs=None, **params):
        procmode = 'sq' if self.mode == 'sq' else None
        caltype = "processed_arc"
        self.getCalibration(adinputs, caltype=caltype,
                            procmode=procmode,
                            refresh=params["refresh"],
                            howmany=params["howmany"])
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedSlit(self, adinputs=None, **params):
        procmode = 'sq' if self.mode == 'sq' else None
        caltype = "processed_slit"
        self.getCalibration(adinputs, caltype=caltype, procmode=procmode,
                            refresh=params["refresh"])
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedSlitFlat(self, adinputs=None, **params):
        procmode = 'sq' if self.mode == 'sq' else None
        caltype = "processed_slitflat"
        self.getCalibration(adinputs, caltype=caltype, procmode=procmode,
                            refresh=params["refresh"])
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    # =========================== STORE PRIMITIVES =================================
    def storeProcessedSlit(self, adinputs=None, **params):
        caltype = 'processed_slit'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self._markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITIM")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedSlitFlat(self, adinputs=None, **params):
        caltype = 'processed_slitflat'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self._markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITFL")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs
