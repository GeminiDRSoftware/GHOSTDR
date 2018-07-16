# pytest suite
"""
Tests for primitives_ghost_bundle.

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


class TestGhostBundle:
    """
    Suite of tests for the functions in the primitives_ghost_bundle module
    """

    def test_splitBundle(self):
        """
        Checks to make:

        - Count number of extracted files
        - Ensure function itself returns empty list
        """
        pass