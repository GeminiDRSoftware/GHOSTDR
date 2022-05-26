#
#                                                                  gemini_python
#
#                                                      primitives_ghost_spect.py
# ------------------------------------------------------------------------------
import sys
import matplotlib.pylab as plt
import astrodata
import ghost_instruments
from ghostdr.ghost.lookups import polyfit_dict
from ghostdr.ghost.polyfit import SlitView

for fname in sys.argv[1:]:
    ad = astrodata.open(fname)
    res_mode = ad.res_mode()
    try:
        slitv_fn = polyfit_dict.get_polyfit_filename(None, 'slitv', res_mode, ad.ut_date(),
                                                     ad.filename, 'slitvmod')
        slitvpars = astrodata.open(slitv_fn)
        print(f"Using slitvmod {slitv_fn}")
    except IOError:
        sys.exit(1)

    sview = SlitView(ad[0].data, None,
                     slitvpars.TABLE[0], mode=res_mode,
                     microns_pix = 4.54 * 180 / 50 * ad.detector_x_bin())
    rcutout = sview.cutout('red')
    rprofile = sview.slit_profile('red')
    bcutout = sview.cutout('blue')
    bprofile = sview.slit_profile('blue')
    f, axarr = plt.subplots(2,2)
    f.suptitle(fname)
    axarr[0,0].imshow(rcutout, origin='lower')
    axarr[0,0].set_title("Red")
    axarr[0,1].imshow(bcutout, origin='lower')
    axarr[0,1].set_title("Blue")
    axarr[1,0].plot(rprofile)
    axarr[1,0].set_title("Red")
    axarr[1,1].plot(bprofile)
    axarr[1,1].set_title("Blue")
plt.show()
