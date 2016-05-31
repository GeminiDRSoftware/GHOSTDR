from GHOST_Keywords import GHOST_KeyDict

from astrodata_Gemini.ADCONFIG_Gemini.descriptors.GEMINI_Descriptors import GEMINI_DescriptorCalc

class GHOST_DescriptorCalc(GEMINI_DescriptorCalc):
    _update_stdkey_dict = GHOST_KeyDict

    def __init__(self):
        GEMINI_DescriptorCalc.__init__(self)


