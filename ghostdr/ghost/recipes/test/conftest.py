import os
import py

FULL_REDUCTION_TMPDIR = 'ghost_fullreduce'

FULL_REDUCTION_SPACE_REQD = 5. * 1024.  # MB

THIS_DIR = os.path.dirname(__file__)

import pytest
import ctypes
import platform
import sys
import glob
import shutil
import subprocess

from recipe_system import __version__
from recipe_system.utils.reduce_utils import buildParser
from recipe_system.utils.reduce_utils import normalize_args
from recipe_system.cal_service import set_calservice
from recipe_system.cal_service import CalibrationService
from recipe_system.config import globalConf


def get_free_space_mb(dirname):
    """Return folder/drive free space (in megabytes)."""
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(dirname),
                                                   None, None,
                                                   ctypes.pointer(free_bytes))
        return free_bytes.value / 1024. / 1024.
    else:
        st = os.statvfs(dirname)
        return st.f_bavail * st.f_frsize / 1024. / 1024.


@pytest.fixture(scope='session')
def get_or_create_tmpdir(tmpdir_factory):
    basetmp = tmpdir_factory.getbasetemp()

    # This test suite requires a minimum amount of available disk space.
    #  Will raise a RuntimeError if this isn't the case.
    if get_free_space_mb(os.path.join(
            basetmp.dirname, basetmp.basename)) < FULL_REDUCTION_SPACE_REQD:
        raise RuntimeError('You have insufficient free disk space to run '
                           'the full reduction test suite.')

    try:
        os.chdir(os.path.join(basetmp.dirname, basetmp.basename,
                              FULL_REDUCTION_TMPDIR))
        tmpsubdir = py.path.local(os.getcwd())
        print('tmpsubdir is {}'.format(tmpsubdir))
    except OSError:
        tmpsubdir = tmpdir_factory.mktemp(FULL_REDUCTION_TMPDIR,
                                          numbered=False)
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))

    # OLD WAY
    # Blank the calibrations manager
    # print('Blanking calibrations manager')
    # subprocess.check_call(['caldb', 'init',  '-v', '-w', ])

    # NEW WAY
    # Set up the calibration system with appropriate arguments
    os.mkdir('dbdir')
    parser = buildParser(__version__)
    local_db_dir = u'{}'.format(os.path.join(tmpsubdir.dirname,
                                             tmpsubdir.basename, 'dbdir/'), )
    args = parser.parse_args(args=[
        '--local_db_dir',
        local_db_dir,
    ])
    args = normalize_args(args)

    # set_calservice(args)
    #
    # # Get a LocalManager and instantiate
    # # import pdb; pdb.set_trace()
    # lm = LocalManager(globalConf['calibs'].database_dir)
    # lm.init_database(wipe=True)

    cs = CalibrationService()
    cs.config(db_dir=local_db_dir)
    cs.init(wipe=True)

    # TODO If necessary, populate calibration system from testdata/calibs
    # Ideally, however, the database would be populated 'as we go' during the
    # tests

    yield tmpsubdir, cs

    # Teardown code - clear out the tmpdir, except the log files
    for _ in glob.glob(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            '*.fits'),
    ):
        os.remove(_)

    # Remove the detritus of the calibrations system
    try:
        shutil.rmtree(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            'calibrations'))
    except OSError:
        pass
    try:
        shutil.rmtree(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            'dbdir'))
    except OSError:
        pass
