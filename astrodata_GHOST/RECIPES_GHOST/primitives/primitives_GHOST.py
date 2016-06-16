from astrodata import AstroData
from astrodata.utils import logutils
from astrodata.utils import Errors
from gempy.gemini import gemini_tools as gt
from astrodata_GHOST.ADCONFIG_GHOST.lookups import timestamp_keywords as ghost_stamps

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
        self.timestamp_keys.update(ghost_stamps.timestamp_keys)
        return rc

    def validateData(self, rc):
        """
        This primitive is used to validate GHOST data, specifically.
        """

        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "validateData", "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["validateData"]

        # Initialize the list of output AstroData objects
        adoutput_list = []

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():

            # Check whether the validateData primitive has been run previously
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by validateData"
                            % ad.filename)

                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # Validate the input AstroData object by ensuring that it has
            # 4 science extensions (NOTE: this might not be correct for
            # spectroscopy)
            num_ext = ad.count_exts("SCI")
            if num_ext != 4:
                raise Errors.Error("The number of extensions in %s do "
                                   "match with the number of extensions "
                                   "expected in raw GHOST data."
                                   % ad.filename)
            else:
                log.fullinfo("The GHOST input file has been validated: %s "
                             "contains %d extensions" % (ad.filename, num_ext))

            # Add the appropriate time stamps to the PHU
            gt.mark_history(adinput=ad, primname=self.myself(), keyword=timestamp_key)

            # Change the filename
            ad.filename = gt.filename_updater(adinput=ad, suffix=rc["suffix"],
                                              strip=True)

            # Append the output AstroData object to the list of output
            # AstroData objects
            adoutput_list.append(ad)

        # Report the list of output AstroData objects to the reduction context
        rc.report_output(adoutput_list)

        yield rc

    def standardizeStructure(self, rc):
        """
        This primitive is used to standardize the structure of GHOST data,
        specifically.

        :param attach_mdf: Set to True to attach an MDF extension to the input
                           AstroData object(s). If an input AstroData object
                           does not have an AstroData type of SPECT, no MDF
                           will be added, regardless of the value of this
                           parameter.
        :type attach_mdf: Python boolean
        :param mdf: The file name, including the full path, of the MDF(s) to
                    attach to the input AstroData object(s). If only one MDF is
                    provided, that MDF will be attached to all input AstroData
                    object(s). If more than one MDF is provided, the number of
                    MDFs must match the number of input AstroData objects. If
                    no MDF is provided, the primitive will attempt to determine
                    an appropriate MDF.
        :type mdf: string or list of strings
        """
        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "standardizeStructure",
                                 "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["standardizeStructure"]

        # Initialize the list of output AstroData objects
        adoutput_list = []

        # Use a flag to determine whether to run addMDF
        attach_mdf = True

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():

            # Check whether the standardizeStructure primitive has been run
            # previously
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by standardizeStructure"
                            % ad.filename)

                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # Attach an MDF to each input AstroData object
            if rc["attach_mdf"] and attach_mdf:

                # Get the mdf parameter from the reduction context
                mdf = rc["mdf"]
                if mdf is not None:
                    rc.run("addMDF(mdf=%s)" % mdf)
                else:
                    rc.run("addMDF")

                # Since addMDF uses all the AstroData inputs from the reduction
                # context, it only needs to be run once in this loop
                attach_mdf = False

            # Add the appropriate time stamps to the PHU
            gt.mark_history(adinput=ad, primname=self.myself(), keyword=timestamp_key)

            # Change the filename
            ad.filename = gt.filename_updater(adinput=ad, suffix=rc["suffix"],
                                              strip=True)

            # Append the output AstroData object to the list of output
            # AstroData objects
            adoutput_list.append(ad)

        # Report the list of output AstroData objects to the reduction context
        rc.report_output(adoutput_list)

        yield rc

    def standardizeHeaders(self, rc):
        """
        This primitive is used to standardize the headers of GHOST data,
        specifically.
        """
        # Instantiate the log
        log = logutils.get_logger(__name__)

        # Log the standard "starting primitive" debug message
        log.debug(gt.log_message("primitive", "standardizeHeaders",
                                 "starting"))

        # Define the keyword to be used for the time stamp for this primitive
        timestamp_key = self.timestamp_keys["standardizeHeaders"]

        # Initialize the list of output AstroData objects
        adoutput_list = []

        # Loop over each input AstroData object in the input list
        for ad in rc.get_inputs_as_astrodata():

            # Check whether the standardizeStructure primitive has been run
            # previously
            if ad.phu_get_key_value(timestamp_key):
                log.warning("No changes will be made to %s, since it has "
                            "already been processed by standardizeHeaders"
                            % ad.filename)

                # Append the input AstroData object to the list of output
                # AstroData objects without further processing
                adoutput_list.append(ad)
                continue

            # Add the appropriate time stamps to the PHU
            gt.mark_history(adinput=ad, primname=self.myself(), keyword=timestamp_key)

            # Change the filename
            ad.filename = gt.filename_updater(adinput=ad, suffix=rc["suffix"],
                                              strip=True)

            # Append the output AstroData object to the list of output
            # AstroData objects
            adoutput_list.append(ad)

        # Report the list of output AstroData objects to the reduction context
        rc.report_output(adoutput_list)

        yield rc

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

