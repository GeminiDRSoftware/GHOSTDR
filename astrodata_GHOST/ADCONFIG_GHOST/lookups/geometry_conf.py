# GHOST instrument geometry configuration parameters for
# gem_mosaic_function() in gemMosaicFunction.py

ubin_red = (
    '../../../astrodata_GHOST/ADCONFIG_GHOST/lookups', 'E2V-CCD-231-C6', None,
    'unbinned')
ubin_blue = (
    '../../../astrodata_GHOST/ADCONFIG_GHOST/lookups', 'E2V-CCD-231-84', None,
    'unbinned')

# gaps: (x_gap, y_gap). Instr, dettype, detector, bin/unbinned
# Gaps between the blocks in the mosaic grid, (0, 0) (col, row) is the lower
# left.
gaps_tile = {
    ubin_red: {(0, 0): (0, 0)},
    ubin_blue: {(0, 0): (0, 0)},
}

# gaps: (x_gap, y_gap) applied when mosaicking corrected blocks. Instr, dettype,
# detector, bin/unbinned
gaps_transform = {
    ubin_red: {(0, 0): (0, 0)},
    ubin_blue: {(0, 0): (0, 0)},
}

blocksize = {  # (y, x)
    ubin_red: (6160, 6144),
    ubin_blue: (4112, 4096),
}

shift = {
    ubin_red: [(0., 0.)],
    ubin_blue: [(0., 0.)],
}

rotation = {  # unbinned, units are Degrees.
    ubin_red: (0.,),
    ubin_blue: (0.,),
}

magnification = {  # (x_mag, y_mag)
     ubin_red: [(1., 1.)],
     ubin_blue: [(1., 1.)],
}

mosaic_grid = {
    ubin_red: (1, 1),
    ubin_blue: (1, 1),
}

interpolator = {  # Values could be 'linear', 'poly3', 'poly5', 'spline3'
    'SCI':  'linear', 'DQ':  'linear', 'VAR':  'linear', 'CR':  'linear',
}

ref_block = {
    'ref_block':  (0, 0),  # (0-based). Reference detector position
                           # (x, y) in the mosaic grid.
}
