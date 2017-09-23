#
#                                                                  gemini_python
#
#                                                    primitives_calibdb_ghost.py
# ------------------------------------------------------------------------------
from gempy.gemini import gemini_tools as gt

from geminidr.core import CalibDB
from .parameters_calibdb_ghost import ParametersCalibDBGHOST

from recipe_system.utils.decorators import parameter_override
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
        self.parameters = ParametersCalibDBGHOST

    def getProcessedSlit(self, adinputs=None, **params):
        caltype = "processed_slit"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedSlitBias(self, adinputs=None, **params):
        caltype = "processed_slitbias"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedSlitDark(self, adinputs=None, **params):
        caltype = "processed_slitdark"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedSlitFlat(self, adinputs=None, **params):
        caltype = "processed_slitflat"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedWavefit(self, adinputs=None, **params):
        caltype = "processed_wavefit"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    def getProcessedXmod(self, adinputs=None, **params):
        caltype = "processed_xmod"
        self.getCalibration(adinputs, caltype=caltype)
        self._assert_calibrations(adinputs, caltype)
        return adinputs

    # =========================== STORE PRIMITIVES =================================
    def storeProcessedSlit(self, adinputs=None, **params):
        caltype = 'processed_slit'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITIM")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedSlitBias(self, adinputs=None, **params):
        caltype = 'processed_slitbias'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITBI")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedSlitDark(self, adinputs=None, **params):
        caltype = 'processed_slitdark'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITDA")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedSlitFlat(self, adinputs=None, **params):
        caltype = 'processed_slitflat'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRSLITFL")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedWavefit(self, adinputs=None, **params):
        caltype = 'processed_wavefit'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRWAVLFT")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs

    def storeProcessedXmod(self, adinputs=None, **params):
        caltype = 'processed_xmod'
        sfx = params["suffix"]
        self.log.debug(gt.log_message("primitive", self.myself(), "starting"))
        adinputs = self.markAsCalibration(adinputs, suffix=sfx,
                                    primname=self.myself(), keyword="PRPOLYFT")
        self.storeCalibration(adinputs, caltype=caltype)
        return adinputs
