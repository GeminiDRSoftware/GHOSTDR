"""
Test suite for GHOSTDR tags
"""

import pytest
import astrodata
import ghost_instruments
import os

THIS_DIR = os.path.dirname(__file__)

from lut_tags import fixture_data as tags_fixture_data

# ---
# REGRESSION TESTING
# ---


class FixtureIterator(object):
    """
    Iterate over all files in a directory, returning each attached tagset.
    """
    def __init__(self, data_dict):
        """
        Parameters
        ----------
        data_dict : dict
            A dictionary, of the form::

                {
                    ('GHOST', 'filename.fits'): ['tag', 'set', 'in', 'full', ],
                    ...
                }

        This dictionary is imported from :any:`test.lut_tags`.
        """
        self._data = data_dict

    def __iter__(self):
        for key in sorted(self._data.keys()):
            (instr, filename) = key
            # ad = astrodata.open(os.path.join(THIS_DIR,
            #                                  'test_data', instr, filename
            #                                  ))
            ad = astrodata.open(os.path.join(
                # '/Users/marc/Documents/ghost/testdata-181121',
                THIS_DIR,
                'testdata',
                filename
            ))
            yield filename, ad, set(self._data[key])


@pytest.mark.parametrize("fn,ad,tag_set", FixtureIterator(tags_fixture_data))
def test_descriptor(fn, ad, tag_set):
    """
    Ensure the tag set returned from each test file is as expected.
    """
    assert ad.tags == tag_set
