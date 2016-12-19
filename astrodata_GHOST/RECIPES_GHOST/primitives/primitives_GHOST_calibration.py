from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors

from gempy.gemini import gemini_tools as gt

from astrodata_Gemini.RECIPES_Gemini.primitives.primitives_calibration import \
    CalibrationPrimitives

class GHOST_CalibrationPrimitives(CalibrationPrimitives):
    """
    Calibration fetch and store primitive set for GHOST.

    Attributes
    ----------
    astrotype : str
        Set to "GHOST"
    """

    astrotype = "GHOST"

    def getProcessedPolyfit(self, rc):
        # Instantiate the log
        log = logutils.get_logger(__name__)

        caltype = "processed_polyfit"
        source = rc["source"]
        if source == None:
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
                        raise Errors.InputError("Calibration not found for " \
                                                "%s" % ad.filename)
                else:
                    if first:
                        log.stdinfo("getCalibration: Results")
                        first = False
                    log.stdinfo("   %s\n      for %s" % (cal.filename,
                                                         ad.filename))
            else:
                if "qa" not in rc.context:
                    raise Errors.InputError("Calibration not found for %s" %
                                            ad.filename)

        yield rc

    def storeProcessedPolyfit(self, rc):
        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "storeProcessedPolyfit",
                                 "starting"))

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():
            # Updating the file name with the suffix for this primitive and
            # then report the new file to the reduction context
            ad.filename = gt.filename_updater(adinput=ad, suffix=rc["suffix"],
                                              strip=False)

            # Adding a PROCPOLY time stamp to the PHU
            gt.mark_history(adinput=ad, primname=self.myself(),
                            keyword="PROCPOLY")

            # Refresh the AD types to reflect new processed status
            ad.refresh_types()

        # Upload polyfit(s) to cal system
        rc.run("storeCalibration")

        yield rc
