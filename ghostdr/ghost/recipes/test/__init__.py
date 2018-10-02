"""
This is the 'all-up' test set for GHOSTDR.

These tests will effectively do a full data reduction of a simulated test
data set, encompassing all arms and resolution modes. The standard recipes
have been broken apart into smaller recipes, so the output of each can be
tested.

In order the force tests to be executed in the correct order, each module
has a numeric prefix (as pytest will execute the test modules in alphanumeric
order):

+--------------+--------+
| Obs. type    | Prefix |
+--------------+--------+
| Bundles      | 00n_   |
+--------------+--------+
| Slit bias    | 01n_   |
+--------------+--------+
| Slit dark    | 02n_   |
+--------------+--------+
| Slit flat    | 03n_   |
+--------------+--------+
| Slit arc     | 04n_   |
+--------------+--------+
| Slit         | 05n_   |
+--------------+--------+
| Bias         | 11n_   |
+--------------+--------+
| Dark         | 12n_   |
+--------------+--------+
| Flat         | 13n_   |
+--------------+--------+
| Arc          | 14n_   |
+--------------+--------+
| Standard     | 15n_   |
+--------------+--------+
| Science      | 16n_   |
+--------------+--------+


"""

import os
import py

FULL_REDUCTION_TMPDIR = 'ghost_fullreduce'

FULL_REDUCTION_SPACE_REQD = 5. * 1024.  # MB

import ctypes
import platform
import sys


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


def get_or_create_tmpdir(tf):
    basetmp = tf.getbasetemp()

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
        tmpsubdir = tf.mktemp(FULL_REDUCTION_TMPDIR,
                              numbered=False)
        os.chdir(os.path.join(tmpsubdir.dirname, tmpsubdir.basename))
    return tmpsubdir
