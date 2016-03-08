from astrodata import AstroData
from gempy.gemini import gemini_tools at gt

from primitives_GEMINI import GEMINIPrimitives

class GHOSTPrimitives(GEMINIPrimitives):
    """
    Class containing all the GHOST primitives.  It inherits all the primitives
    from the GEMINIPrimitives class (which itself inherits a series of 
    primitives from RECIPE_Gemini/primitives.)
    """

    astrotype = "GHOST"

    def init(self, rc):
        GEMINIPrimitives.init(self, rc)
        return rc

    def myScienceStep(self, rc):
        """
        Describe primitive.
        """

        # Setup the logger
        log = logutils.get_logger(__name__)
        log.debug(gt.log_message("primitive", "myScienceStep", "starting"))

        # Define the keyword to be used for the timestamp for this primitive
        timestamp_key = self.timestamp_keys['MYSCISTP']

        # Initialize the list of output AstroData objects
        adoutput_list = []

        # Loop over each input AstroData object in the input list
        for ad in rc.get_intputs_as_astrodata():
        
             # Skip if already run through this primitive
             # Just put ad in the output list without processing.
             if ad.phu_get_key_value(timestamp_key):
                 log.warning("No changes will be made to %s, since it has "
                            "already been processed by myScienceStep"
                            % ad.filename)
                 adoutput_list.append(ad)
                 continue

             ########
             ########
             #  The code that does the real work goes here
             #  Can be call to functions, call to other primitives
             ########
             ########

             # Add timestamp to the PHU
             gt.mark_history(adinput=ad, keyword=timestamp_key)

             # Assign suffix to output ad object filename
             ad.filename = gt.filename_updater(adinput=ad, suffix=rc['suffix'],
                                               strip=True)

             # Append the processed file to the output list
             adoutput_list.append(ad)

        # Report the outputs to the Reduction Context, rc.
        rc.report_output(adoutput_list)

        # Return control to the Recipe System
        yield rc

