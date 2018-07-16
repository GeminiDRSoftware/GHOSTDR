# pytest suite
"""
Tests for primitives_ghost_slit.

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
    Suite of tests for the functions in the primitives_ghost_slit module
    """

    def test_CRCorrect(self):
        """
        Checks to make:

        - Check that all simulated CR are removed from test data
        - Check shape of output data matches shape of input
        - Checks on the data median/std dev?
        """
        pass


    def test_processSlits(self):
        """
        Checks to make:

        - Ensure the slit viewer bundle ends up with
            a) A mean exposure epoch
            b) The correct mean exposure epoch
        """


    def test_stackFrames(self):
        """
        Checks to make:

        - Only one file comes out
        - Dimensions of the output image match those of the input image
        - Attempting to stack images of different shapes should fail
        """


    def test__mad(self):
        """
        Checks to make:

        - Pass in some known data, check the MAD is computed correctly
        - Check across axes as well
        """


    def test__total_obj_flux(self):
        """
        Checks to make

        - Compare against already-known total flux?
        """