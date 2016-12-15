from GHOST_Keywords import GHOST_KeyDict

from datetime import datetime
from time import strptime

from astrodata.utils import Errors
from astrodata.interface.Descriptors import DescriptorValue

from gempy.gemini import gemini_metadata_utils as gmu

from astrodata_Gemini.ADCONFIG_Gemini.descriptors.GEMINI_Descriptors import GEMINI_DescriptorCalc

class GHOST_DescriptorCalc(GEMINI_DescriptorCalc):
    _update_stdkey_dict = GHOST_KeyDict

    def __init__(self):
        GEMINI_DescriptorCalc.__init__(self)

    def arm(self, dataset, **args):
        arm_col = None

        data_types = dataset.types
        if "GHOST_BLUE" in data_types:
            arm_col = 'blue'
        elif "GHOST_RED" in data_types:
            arm_col = 'red'

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(arm_col, name="arm", ad=dataset)

        return ret_dv

    def gain(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_gain_dict = {}

        # If the data have been prepared, take the gain value directly from the
        # appropriate keyword. At some point, a check for the STDHDRSI header
        # keyword should be added, since the function that overwrites the gain
        # keyword also writes the STDHDRSI keyword.
        if "PREPARED" in dataset.types:

            # Determine the gain keyword from the global keyword dictionary
            keyword = self.get_descriptor_key("key_gain")

            # Get the value of the gain keyword from the header of each pixel
            # data extension as a dictionary where the key of the dictionary is
            # an ("*", EXTVER) tuple
            gain_dict = gmu.get_key_value_dict(adinput=dataset,
                                               keyword=keyword)
            if gain_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info

            ret_gain_dict = gain_dict

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_gain_dict, name="gain", ad=dataset)

        return ret_dv

    def group_id(self, dataset, **args):
        # Adaptation of other instruments' (e.g. NIRI, F2 etc) group_id
        # descriptor to GHOST

        # Descriptors applicable to all data types
        unique_id_descriptor_list_all = ["arm"]

        # List to format descriptor calss using 'pretty=True' parameter
        call_pretty_version_list = []

        # Descriptors to be returned as ordered_list using descriptor
        # 'as_list' method
        convert_to_list_list = []

        # Additional descriptors for each frame type
        dark_id = ["exposure_time", "coadds"]
        flat_twilight_id = ["observation_id"]
        science_id = ["observation_id"]

        # Additional descriptors required for spectra (i.e. actual observations)
        required_spectra_descriptors = []
        if "GHOST_SPECT" in dataset.types:
            unique_id_descriptor_list_all.extend(required_spectra_descriptors)

        # Choose/update the list of descriptors to be used to generate group_id
        # based on the specific image type
        data_types = dataset.types
        id_descriptor_list = []
        if "GHOST_BIAS" in data_types:
            id_descriptor_list = []
        elif "GHOST_DARK" in data_types:
            id_descriptor_list = dark_id
        else:
            id_descriptor_list = science_id

        # Add in all of the common descriptors required
        id_descriptor_list.extend(unique_id_descriptor_list_all)

        # Form the group_id
        descriptor_object_string_list = []
        for descriptor in id_descriptor_list:
            # Prepare the descriptor call
            if descriptor in call_pretty_version_list:
                end_parameter = "(pretty=True)"
            else:
                end_parameter = "()"
            descriptor_call = ''.join([descriptor, end_parameter])

            # Call the descriptor
            exec ("descriptor_object = dataset.{0}".format(descriptor_call))

            # Check for a returned descriptor value object with a None value
            if descriptor_object.is_none():
                # The descriptor functions return None if a value cannot be
                # found and stores the exception info. Re-raise the exception.
                # It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info

            # In some cases require the information as a list
            if descriptor in convert_to_list_list:
                descriptor_object = descriptor_object.as_list()

            # Convert DV value to a string and store
            descriptor_object_string_list.append(str(descriptor_object))

        # Create the final group_id string
        ret_group_id = '_'.join(descriptor_object_string_list)

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_group_id, name="group_id", ad=dataset)

        return ret_dv

    def res_mode(self, dataset, **args):
        # Get the GHOST resolution mode of this dataset
        # Valid outputs are 'std' and 'high' (or None if not applicable)

        mode_val = None

        # Determine the mode keyword from the global keyword dict
        keyword = self.get_descriptor_key("key_res_mode")

        mode_val = dataset.phu_get_key_value(keyword)

        if mode_val == 'HI_ONLY':
            return_mode_val = 'high'
        elif mode_val == 'LO_ONLY':
            return_mode_val = 'std'
        else:
            return_mode_val = None

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(return_mode_val, name="res_mode", ad=dataset)

        return ret_dv

    def read_noise(self, dataset, **args):
        # Since this descriptor function accesses keywords in the headers of
        # the pixel data extensions, always construct a dictionary where the
        # key of the dictionary is an (EXTNAME, EXTVER) tuple
        ret_read_noise_dict = {}

        # If the data have been prepared, take the read noise value directly
        # from the appropriate keyword. At some point, a check for the STDHDRSI
        # header keyword should be added, since the function that overwrites
        # the read noise keyword also writes the STDHDRSI keyword.
        if "PREPARED" in dataset.types:

            # Determine the read noise keyword from the global keyword
            # dictionary
            keyword = self.get_descriptor_key("key_read_noise")

            # Get the value of the read noise keyword from the header of each
            # pixel data extension as a dictionary where the key of the
            # dictionary is an ("*", EXTVER) tuple
            read_noise_dict = gmu.get_key_value_dict(adinput=dataset,
                                                     keyword=keyword)
            if read_noise_dict is None:
                # The get_key_value_dict() function returns None if a value
                # cannot be found and stores the exception info. Re-raise the
                # exception. It will be dealt with by the CalculatorInterface.
                if hasattr(dataset, "exception_info"):
                    raise dataset.exception_info

            for ext_name_ver, raw_read_noise in read_noise_dict.iteritems():
                if raw_read_noise is None:
                    read_noise = None
                else:
                    read_noise = float(raw_read_noise)

                # Update the dictionary with the read noise value
                ret_read_noise_dict.update({ext_name_ver: read_noise})

        # Instantiate the return DescriptorValue (DV) object
        ret_dv = DescriptorValue(ret_read_noise_dict, name="read_noise",
                                 ad=dataset)
        return ret_dv
