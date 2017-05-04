import inspect
from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors

from gempy.gemini import gemini_tools as gt

from primitives_calibration import CalibrationPrimitives

class GHOST_CalibrationPrimitives(CalibrationPrimitives):
    """
    Calibration fetch and store primitive set for GHOST.

    Attributes
    ----------
    astrotype : str
        Set to "GHOST"
    """

    astrotype = "GHOST"

    def getProcessedSlit(self, rc):
        self._getProcessed(rc)
        yield rc

    def getProcessedSlitBias(self, rc):
        self._getProcessed(rc)
        yield rc

    def getProcessedSlitDark(self, rc):
        self._getProcessed(rc)
        yield rc

    def getProcessedSlitFlat(self, rc):
        self._getProcessed(rc)
        yield rc

    def getProcessedXmod(self, rc):
        self._getProcessed(rc)
        yield rc

    def getProcessedWavefit(self, rc):
        self._getProcessed(rc)
        yield rc

    def _getProcessed(self, rc):
        """
        Helper function - prevents code duplication for the variety of
        getProcessed* functions

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Returns
        -------
        Nil. The calling function is responsible for returning the
        ReductionContext dict.

        """
        # Instantiate the log
        log = logutils.get_logger(__name__)

        caltype = inspect.currentframe().f_back.f_code.co_name
        caltype = caltype.split('getProcessed')[-1].lower()
        caltype = 'processed_' + caltype

        source = rc["source"]
        if source is None:
            rc.run("getCalibration(caltype=%s)" % caltype)
        else:
            rc.run("getCalibration(caltype=%s, source=%s)" % (caltype, source))

        # List calibrations found
        first = True
        for ad in rc.get_inputs_as_astrodata():
            calurl = rc.get_cal(ad, caltype)  # get from cache
            if calurl:
                cal = AstroData(calurl)
                if cal.filename is None:
                    if "qa" not in rc.context:
                        raise Errors.InputError("Calibration not found for "
                                                "%s" % ad.filename)
                else:
                    if first:
                        log.stdinfo("getCalibration: Results")
                        first = False
                    log.stdinfo(
                        "   %s\n      for %s" % (cal.filename, ad.filename))
            else:
                if "qa" not in rc.context:
                    raise Errors.InputError(
                        "Calibration not found for %s" % ad.filename)

    def storeProcessedSlit(self, rc):
        self._storeProcessed(rc, 'PRSLITIM')
        yield rc

    def storeProcessedSlitBias(self, rc):
        self._storeProcessed(rc, 'PRSLITBI')
        yield rc

    def storeProcessedSlitDark(self, rc):
        self._storeProcessed(rc, 'PRSLITDA')
        yield rc

    def storeProcessedSlitFlat(self, rc):
        self._storeProcessed(rc, 'PRSLITFL')
        yield rc

    def storeProcessedXmod(self, rc):
        self._storeProcessed(rc, 'PRPOLYFT')
        yield rc

    def storeProcessedWavefit(self, rc):
        self._storeProcessed(rc, 'PRWAVLFT')
        yield rc

    def _storeProcessed(self, rc, key):
        """
        Helper function - prevents code duplication for the variety of
        storeProcessed* functions

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Returns
        -------
        Nil. The calling function is responsible for returning the
        ReductionContext dict.

        """

        # can't use self.myself(): it returns "_storeProcessed" here; instead,
        # we want it to return the caller of _storeProcessed (_storeProcessed
        # is just a delegate)
        me = inspect.currentframe().f_back.f_code.co_name

        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", me, "starting"))

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():
            # Updating the file name with the suffix for this primitive and
            # then report the new file to the reduction context
            ad.filename = gt.filename_updater(adinput=ad, suffix=rc["suffix"],
                                              strip=False)

            # Adding a time stamp to the PHU
            gt.mark_history(adinput=ad, primname=me, keyword=key)

            # Refresh the AD types to reflect new processed status
            ad.refresh_types()

        # Upload the calibration to cal system
        rc.run("storeCalibration")
