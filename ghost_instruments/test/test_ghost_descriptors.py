import pytest
import astrodata
import gemini_instruments
import ghost_instruments
import os

THIS_DIR = os.path.dirname(__file__)

from lut_descriptors import fixture_data as descriptors_fixture_data

# ---
# REGRESSION TESTING
# ---


class FixtureIterator(object):
    """
    Helper class to iterate over all files in a directory, returning an
    attached descriptor and its value in turn.
    """
    def __init__(self, data_dict):
        self._data = data_dict

    def __iter__(self):
        for (instr, filename) in sorted(self._data.keys()):
            # ad = astrodata.open(os.path.join(THIS_DIR,
            #                                  'test_data', instr, filename
            #                                  ))
            ad = astrodata.open(os.path.join(
                # '/Users/marc/Documents/ghost/testdata-181121',
                THIS_DIR,
                'testdata',
                filename
            ))
            for desc, value in self._data[(instr, filename)]:
                yield filename, ad, desc, value


@pytest.mark.parametrize("fn, ad, descriptor, value",
                         FixtureIterator(descriptors_fixture_data))
def test_descriptor(fn, ad, descriptor, value):
    """
    Ensure that the values returned by AstroData descriptors are as expected
    """
    method = getattr(ad, descriptor)
    mvalue = method()
    if float in (type(value), type(mvalue)):
        assert abs(mvalue - value) < 0.0001
    else:
        assert mvalue == value
