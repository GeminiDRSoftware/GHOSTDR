# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost.py file, in alphabetical order.

from geminidr.gemini.parameters_gemini import ParametersGemini
from geminidr.core.parameters_ccd import ParametersCCD
from .parameters_calibdb_ghost import ParametersCalibDBGHOST

class ParametersGHOST(ParametersGemini, ParametersCCD, ParametersCalibDBGHOST):
    pass
