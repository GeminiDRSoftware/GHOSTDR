#
# This script is used to help verify reqirement 5625.
# It should be run on bundled FITS files as downloaded
# from the archive.
#

import sys
from astropy.io import fits

#
# The sets of keywords expected in the various extensions
#
phu = {
    'INSTRUME',
    'OBJECT',
    'OBSTYPE',
    'OBSCLASS',
    'GEMPRGID',
    'OBSID',
    'DATALAB',
    'RA',
    'DEC',
    'EQUINOX',
    'CRPA',
    'CRFOLLOW',
    'ELEVATIO',
    'AZIMUTH',
    'HA',
    'PA',
    'IAA',
    'AIRMASS',
    'AMSTART',
    'AMEND',
    'HUMIDITY',
    'TAMBIENT',
    'PRESSUR2',
    'SIMPLE',
    'BITPIX',
    'NAXIS',
    'EXTEND',
    'NEXTEND',
    'NREDEXP',
    'NBLUEEXP',
    'NSLITEXP',
}

sciphu = {
    'DATE-OBS',
    'UTSTART',
    'UTEND',
    'DARKTIME',
    'DETSIZE',
    'EXPUTST',
    'EXPUTEND',
    'DSPTIMBN',
    'DSPTIMBV'
}

sci = {
    'XTENSION',
    'BITPIX',
    'NAXIS',
    'NAXIS1',
    'NAXIS2',
    'PCOUNT',
    'GCOUNT',
    'BZERO',
    'CAMERA',
    'EXPID',
    'DETSEC',
    'CCDSIZE',
    'CCDNAME',
    'AMPNAME',
    'AMPSIZE',
    'CCDSEC',
    'BIASSEC',
    'TRIMSEC',
    'CCDSUM',
}

slitv = {
    'XTENSION',
    'BITPIX',
    'NAXIS',
    'NAXIS1',
    'NAXIS2',
    'PCOUNT',
    'GCOUNT',
    'BZERO',
    'CAMERA',
    'EXPID',
    'DETSIZE',
    'CCDSIZE',
    'AMPSIZE',
    'DETSEC',
    'CCDSEC',
    'CCDNAME',
    'CCDSUM',
    'EXPUTST',
    'EXPUTEND'
}

# Open every file, iterate over all extensions, match up expected keywords with actual keywords
for fn in sys.argv[1:]:
    with fits.open(fn) as hdus:
        print("File: ", fn)
        hdus_sets = [set(h.header) for h in hdus]
        if not phu.issubset(hdus_sets[0]):
            print("Missing keywords in PHU:", phu.difference(hdus_sets[0]))
        for i, (hdu, s) in enumerate(zip(hdus[1:], hdus_sets[1:])):
            if not 'CAMERA' in hdu.header:
                print("Missing CAMERA keyword in extension", i)
            elif hdu.header['CAMERA'] == 'RED' or hdu.header['CAMERA'] == 'BLUE':
                if hdu.header['NAXIS'] == 0:
                    if not sciphu.issubset(s):
                        print("Missing SCI PHU keywords in extension", i, ":", sciphu.difference(s))
                    pass
                else:
                    if not sci.issubset(s):
                        print("Missing SCI keywords in extension", i, ":", sci.difference(s))
            elif hdu.header['CAMERA'] == 'SLITV':
                if not slitv.issubset(s):
                    print("Missing SLITV keywords in extension", i, ":", slitv.difference(s))
