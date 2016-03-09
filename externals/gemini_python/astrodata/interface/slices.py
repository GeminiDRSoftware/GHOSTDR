from pyfits import ImageHDU
from astrodata import AstroData

def pixel_exts(ad):
    """
    Function receives an AstroData instance, searches HDUs for
    Image extensions (a.k.a., pixel extensions).

    :parameter ad: Astrodata instance
    :type ad:  <AstroData>

    :return: list of Image extensions
    :rtype: <list>

    """
    try:
        assert isinstance(ad, AstroData)
    except AssertionError:
        raise TypeError("argument must be an AstroData object")

    return [ext for ext in ad.hdulist if isinstance(ext, ImageHDU)]
