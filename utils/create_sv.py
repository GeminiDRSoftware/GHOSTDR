from astropy.table import Table
from astropy.io import fits
import astrodata

'''
rota, rotyc, rotxc, yc, xc, extract_hw, skypix0, skypix1, obj0pix0, obj0pix1, obj1pix0, obj1pix1
std,  red,  80.96, 900, 780, 781, 985, 3, 47, 63,   3, 46, 64, 107
std,  blue, 80.96, 900, 780, 771, 946, 3, 47, 63,   3, 46, 64, 107
high, red,  80.96, 900, 780, 770, 857, 2, 82, 106, 11, 81, 4,  9  
high, blue, 80.96, 900, 780, 760, 819, 3, 82, 106, 11, 81, 4,  9  
'''

'''
Restructure into two files:
GHOST_1_1_slitv_std_220501.fits and
GHOST_1_1_slitv_high_220501.fits

Each will contain
rota
rotyc
rotxc
center_y_red
center_x_red
center_y_blue
center_x_blue
ext_hw
skypix0
skypix1
obj0pix0
obj0pix1
obj1pix0
obj1pix1
'''

std = {
          'rota':          [0.0,],
          'rotyc':         [0,],
          'rotxc':         [0,],
          #'rota':          [80.96,],
          #'rotyc':         [900,],
          #'rotxc':         [780,],
          'center_y_red':  [77,],
          'center_x_red':  [65,],
          'center_y_blue': [77,],
          'center_x_blue': [156,],
          #'center_y_red':  [781,],
          #'center_x_red':  [985,],
          #'center_y_blue': [771,],
          #'center_x_blue': [946,],
          'ext_hw':        [3,],
          'skypix0':       [47,],
          'skypix1':       [63,],
          'obj0pix0':      [3,],
          'obj0pix1':      [46,],
          'obj1pix0':      [64,],
          'obj1pix1':      [107,],
      }

high = {
          'rota':          [0.0,],
          'rotyc':         [0,],
          'rotxc':         [0,],
          #'rota':          [80.96,],
          #'rotyc':         [900,],
          #'rotxc':         [780,],
            #'red': [78, 95],
            #'blue': [78, 4]
          'center_y_red':  [78,],
          'center_x_red':  [95,],
          'center_y_blue': [78],
          'center_x_blue': [4],
          #'center_y_red':  [770,],
          #'center_x_red':  [857,],
          #'center_y_blue': [760],
          #'center_x_blue': [819],
          'ext_hw':        [2,],
          'skypix0':       [82,],
          'skypix1':       [106,],
          'obj0pix0':      [11,],
          'obj0pix1':      [81,],
          'obj1pix0':      [4,],
          'obj1pix1':      [9,],
      }

t_std = Table(std)
phu = fits.PrimaryHDU()
ad = astrodata.create(phu)
ad.append(t_std, name='TABLE')
ad.write('std.fits', overwrite=True)

t_high = Table(high)
phu = fits.PrimaryHDU()
ad = astrodata.create(phu)
ad.append(t_high, name='TABLE')
ad.write('high.fits', overwrite=True)
