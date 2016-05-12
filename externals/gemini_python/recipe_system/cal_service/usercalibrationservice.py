#
#                                                                  gemini_python
#
#
#                                                      usercalibrationservice.py
# ------------------------------------------------------------------------------
# $Id: usercalibrationservice.py 5157 2015-02-20 19:21:11Z kanderson $
# ------------------------------------------------------------------------------
__version__      = '$Revision: 5157 $'[11:-2]  # Changed by swapper, 22 May 2014
__version_date__ = '$Date: 2015-02-21 06:21:11 +1100 (Sat, 21 Feb 2015) $'[7:-3]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
class UserCalibrationService(object):
    
    def __init__(self):
        self.user_cal_dict = {}
    
    def add_calibration(self, caltype=None, cal_file=None):
        self.user_cal_dict.update({caltype:cal_file})
        return
        
    def get_calibration(self, caltype=None):
        return self.user_cal_dict.get(caltype)

user_cal_service = UserCalibrationService()
