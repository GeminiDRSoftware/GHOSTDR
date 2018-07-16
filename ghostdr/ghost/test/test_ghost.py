# pytest suite
"""
Tests for primitives_ghost.

This is a suite of tests to be run with pytest.

To run:
    1) Set the environment variable GEMPYTHON_TESTDATA to the path that
       contains the directories with the test data.
       Eg. /net/chara/data2/pub/gempython_testdata/
    2) From the ??? (location): pytest -v --capture=no
"""
import os
import numpy as np
import astrodata
import gemini_instruments
from gempy.utils import logutils

# from geminidr.core.test import ad_compare

TESTDATAPATH = os.getenv('GEMPYTHON_TESTDATA', '.')
logfilename = 'test_standardize.log'


class TestGhost:
    """
    Suite of tests for the functions in the ghost_primitives module
    """

    def test__rebing_ghost_ad(self):
        """
        Checks to make:

        - Re-binned data arrays are the correct shape;
        - Correct keyword headers have been updated;

        Loop over each valid binning mode
        """
        pass