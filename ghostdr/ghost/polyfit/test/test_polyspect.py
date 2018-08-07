from __future__ import division, print_function
import pytest
import ghostdr.ghost.polyfit as polyfit
import astropy.io.fits as pyfits
import pdb
import numpy as np

# Create a generic instantiation of ghost for generic test purposes.
# gen_ps = polyfit.polyspect.Polyspect()

# Assert if all the correct attributes of the ghost class needed are there

def test_polyspect_init():
    with pytest.raises(TypeError):
        # Must provide all init parameters
        _ = polyfit.polyspect.Polyspect()

@pytest.mark.skip
@pytest.mark.parametrize("attrib, tp", [
    ('m_ref', int),
    ('szx', int),
    ('szy', int),
    ('m_min', int),
    ('m_max', int),
    ('transpose', bool),
])
def test_polyspect_attributes(attrib, tp):
    """
    Test if the ghost class contains all the needed attributes
    in the right format
    """
    assert hasattr(gen_ps, attrib), "GhostArm missing attribute " \
                                       "{}".format(attrib)
    assert isinstance(getattr(gen_ps, attrib),
                      tp), "GhostArm attribute {} has incorrect type " \
                           "(expected {}, got {})".format(
        attrib, tp, type(getattr(gen_ps, attrib)),
    )
