import os
import py

FULL_REDUCTION_TMPDIR = 'ghost_fullreduce'

FULL_REDUCTION_SPACE_REQD = 5. * 1024.  # MB

import pytest
import ctypes
import platform
import sys
import glob
import shutil


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

    yield tmpsubdir

    # Teardown code - clear out the tmpdir, except the log files
    for _ in glob.glob(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            '*.fits'),
    ):
        os.remove(_)
    try:
        shutil.rmtree(os.path.join(
            tmpsubdir.dirname, tmpsubdir.basename,
            'calibrations'))
    except OSError:
        pass