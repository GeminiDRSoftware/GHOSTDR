# This is a dictionary containing the assignments of initial Polyfit models
# to various GHOST frame types

# Keys are of the format
# instrument_detectorXBin_detectorYBin_Arm_ResMode
import os
from datetime import datetime

xmod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_red_std_161120': 'red/std/161120/xmod.fits',
    'GHOST_1_1_red_high_161120': 'red/high/161120/xmod.fits',
    'GHOST_1_1_blue_std_161120': 'blue/std/161120/xmod.fits',
    'GHOST_1_1_blue_high_161120': 'blue/high/161120/xmod.fits',
    'GHOST_1_1_red_std_220501': 'red/std/220501/xmod.fits',
    'GHOST_1_1_red_high_220501': 'red/high/220501/xmod.fits',
    'GHOST_1_1_blue_std_220501': 'blue/std/220501/xmod.fits',
    'GHOST_1_1_blue_high_220501': 'blue/high/220501/xmod.fits',
    #'GHOST_1_1_blue_high_170801': 'blue/high/170801/xmod.fits',
}

wavemod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_red_std_161120': 'red/std/161120/wavemod.fits',
    'GHOST_1_1_red_high_161120': 'red/high/161120/wavemod.fits',
    'GHOST_1_1_blue_std_161120': 'blue/std/161120/wavemod.fits',
    'GHOST_1_1_blue_high_161120': 'blue/high/161120/wavemod.fits',
    'GHOST_1_1_red_std_220501': 'red/std/220501/wavemod.fits',
    'GHOST_1_1_red_high_220501': 'red/high/220501/wavemod.fits',
    'GHOST_1_1_blue_std_220501': 'blue/std/220501/wavemod.fits',
    'GHOST_1_1_blue_high_220501': 'blue/high/220501/wavemod.fits',
}

spatmod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_red_std_161120': 'red/std/161120/spatmod.fits',
    'GHOST_1_1_red_high_161120': 'red/high/161120/spatmod.fits',
    'GHOST_1_1_blue_std_161120': 'blue/std/161120/spatmod.fits',
    'GHOST_1_1_blue_high_161120': 'blue/high/161120/spatmod.fits',
    'GHOST_1_1_red_std_220501': 'red/std/220501/spatmod.fits',
    'GHOST_1_1_red_high_220501': 'red/high/220501/spatmod.fits',
    'GHOST_1_1_blue_std_220501': 'blue/std/220501/spatmod.fits',
    'GHOST_1_1_blue_high_220501': 'blue/high/220501/spatmod.fits',
}

specmod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_red_std_161120': 'red/std/161120/specmod.fits',
    'GHOST_1_1_red_high_161120': 'red/high/161120/specmod.fits',
    'GHOST_1_1_blue_std_161120': 'blue/std/161120/specmod.fits',
    'GHOST_1_1_blue_high_161120': 'blue/high/161120/specmod.fits',
    'GHOST_1_1_red_std_220501': 'red/std/220501/specmod.fits',
    'GHOST_1_1_red_high_220501': 'red/high/220501/specmod.fits',
    'GHOST_1_1_blue_std_220501': 'blue/std/220501/specmod.fits',
    'GHOST_1_1_blue_high_220501': 'blue/high/220501/specmod.fits',
}

rotmod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_red_std_161120': 'red/std/161120/rotmod.fits',
    'GHOST_1_1_red_high_161120': 'red/high/161120/rotmod.fits',
    'GHOST_1_1_blue_std_161120': 'blue/std/161120/rotmod.fits',
    'GHOST_1_1_blue_high_161120': 'blue/high/161120/rotmod.fits',
    'GHOST_1_1_red_std_220501': 'red/std/220501/rotmod.fits',
    'GHOST_1_1_red_high_220501': 'red/high/220501/rotmod.fits',
    'GHOST_1_1_blue_std_220501': 'blue/std/220501/rotmod.fits',
    'GHOST_1_1_blue_high_220501': 'blue/high/220501/rotmod.fits',
}

slitvmod_dict = {
    # 'key': 'relative_dir',
    'GHOST_1_1_slitv_std_161120': 'slitv/std/161120/slitvmod.fits',
    'GHOST_1_1_slitv_high_161120': 'slitv/high/161120/slitvmod.fits',
    'GHOST_1_1_slitv_std_220501': 'slitv/std/220501/slitvmod.fits',
    'GHOST_1_1_slitv_high_220501': 'slitv/high/220501/slitvmod.fits',
}

def get_polyfit_filename(log, arm, mode, date_obs, filename, caltype):
    """
    Gets the filename of the relevant initial polyfit file for this
    input GHOST science image

    This primitive uses the arm, resolution mode and observing epoch
    of the input AstroData object to determine the correct initial
    polyfit model to provide. The model provided matches the arm and
    resolution mode of the data, and is the most recent model generated
    before the observing epoch.

    Parameters
    ----------
    log: The recipe logger
    arm: The GHOST arm (blue, red, or slitv)
    mode: The GHOST mode (std or high)
    date_obs: ad.ut_date()
    filename: ad.filename
    caltype : str
        The initial model type (e.g. ``'rotmod'``, ``'spatmod'``, etc.)
        requested. An :any:`AttributeError` will be raised if the requested
        model type does not exist.

    Returns
    -------
    str/None:
        Filename (including path) of the required polyfit file
    """
    polyfit_dir = os.path.join(os.path.dirname(__file__),
                               'Polyfit')

    # CJS: This is a method that only exists *if* the input is of type
    # GHOST, so no need to check
    key = 'GHOST_1_1_{}_{}'.format(arm, mode)

    try:
        poly_dict = globals()['{}_dict'.format(caltype)]
    except AttributeError:
        if not log is None:
            log.warning("Invalid polyfit calibration type ({}) requested for "
                        "{}".format(caltype, filename))
        return None

    #FIXME: Restrict search to correct res and arm (all not necessarily
    #updated at once!)
    dates_avail = set([k.split('_')[-1] for k in poly_dict.keys()])
    
    # Safe to assume instrument won't be used after 2099...
    dates_avail = [datetime.strptime('20{}'.format(x),
                                        '%Y%m%d').date() for x in dates_avail]
    dates_avail.sort()

    # Determine the latest data that precedes the observing date
    try:
        date_req = [_ for _ in dates_avail if _ <= date_obs][-1]
    except IndexError:
        if not log is None:
            log.warning("No polyfit available for {}".format(filename))
        return None
    key += '_{}'.format(date_req.strftime('%y%m%d'))

    polyfit_file = poly_dict[key]
    # Prepend standard path if the filename doesn't start with '/'
    return polyfit_file if polyfit_file.startswith(os.path.sep) else \
        os.path.join(polyfit_dir, polyfit_file)

