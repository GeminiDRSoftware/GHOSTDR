from astropy.table import Table
from astropy.io import fits
import astrodata

std = {
          # For simulated data
          'rota':          [0.0,],
          'rotyc':         [0,],
          'rotxc':         [0,],
          'center_y_red':  [154,],
          'center_x_red':  [130,],
          'center_y_blue': [154,],
          'center_x_blue': [312,],
          # For lab data 20220501
          #'rota':          [80.96,],
          #'rotyc':         [900,],
          #'rotxc':         [780,],
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
          # For simulated data
          'rota':          [0.0,],
          'rotyc':         [0,],
          'rotxc':         [0,],
          'center_y_red':  [156,],
          'center_x_red':  [190,],
          'center_y_blue': [156],
          'center_x_blue': [8],
          # For lab data 20220501
          #'rota':          [80.96,],
          #'rotyc':         [900,],
          #'rotxc':         [780,],
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
