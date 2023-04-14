import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
from astropy.wcs import WCS
from astropy.io import fits
from sklearn.metrics.pairwise import euclidean_distances

# Extract from GAIA DR3 for selected objects
# TODO this could be downloaded programmatically
objects = {
        'HD 131977': { "_V": "VizieR", "RA_ICRS": "224.37159427348", "DE_ICRS": "-21.42314039626", "Source": "6232511606838403968", "e_RA_ICRS": "0.0539", "e_DE_ICRS": "0.0468", "Plx": "169.8843", "e_Plx": "0.0653", "PM": "2008.680", "pmRA": "1031.472", "e_pmRA": "0.068", "pmDE": "-1723.619", "e_pmDE": "0.055", "RUWE": "1.298", "FG": "134688970.74624", "e_FG": "9.1027e+04", "Gmag": "5.364037", "FBP": "5.6385e+07", "e_FBP": "7.5198e+04", "BPmag": "5.960628", "FRP": "1.1075e+08", "e_FRP": "1.6761e+05", "RPmag": "4.637024", "BP-RP": "1.323604", "RV": "26.75", "e_RV": "0.12", "Vbroad": "6.9112", "GRVSmag": "4.383033", "QSO": "0", "Gal": "0", "NSS": "0", "XPcont": "0", "XPsamp": "0", "RVS": "0", "EpochPh": "0", "EpochRV": "0", "MCMCGSP": "1", "MCMCMSC": "0", "And": "0", "Teff": "4503.0", "logg": "4.5357", "[Fe/H]": "-0.0429", "Dist": "5.8870", "A0": "0.0000", "HIP": "73184", "PS1": " ", "SDSS13": " ", "SKYM2": "103698048", "TYC2": "6180-855-1", "URAT1": " ", "AllWISE": " ", "APASS9": "56337346", "GSC23": "S91P000642", "RAVE5": " ", "2MASS": "14572788-2124526", "RAVE6": " ", "RAJ2000": "224.36666959995", "DEJ2000": "-21.41547922619" },
        'HD60753':   { "_V": "VizieR", "RA_ICRS": "113.36380749102", "DE_ICRS": "-50.58422947157", "Source": "5493730399606077440", "e_RA_ICRS": "0.0331", "e_DE_ICRS": "0.0321", "Plx": "1.4719", "e_Plx": "0.0342", "PM": "6.320", "pmRA": "-2.897", "e_pmRA": "0.041", "pmDE": "5.617", "e_pmDE": "0.038", "RUWE": "0.975", "FG": "40958791.40509", "e_FG": "1.6622e+04", "Gmag": "6.656499", "FBP": "3.1110e+07", "e_FBP": "3.6438e+04", "BPmag": "6.606296", "FRP": "1.6607e+07", "e_FRP": "1.1778e+04", "RPmag": "6.697140", "BP-RP": "-0.090845", "RV": " ", "e_RV": " ", "Vbroad": " ", "GRVSmag": " ", "QSO": "0", "Gal": "0", "NSS": "0", "XPcont": "1", "XPsamp": "1", "RVS": "0", "EpochPh": "0", "EpochRV": "0", "MCMCGSP": "1", "MCMCMSC": "1", "And": "0", "Teff": "14602.5", "logg": "3.3036", "[Fe/H]": "-1.0603", "Dist": "694.4103", "A0": "0.2154", "HIP": "36745", "PS1": " ", "SDSS13": " ", "SKYM2": "470198731", "TYC2": "8141-1632-1", "URAT1": " ", "AllWISE": "J073327.31-503503.2", "APASS9": "42998650", "GSC23": "S602043235", "RAVE5": " ", "2MASS": "07332733-5035032", "RAVE6": " ", "RAJ2000": "113.36382777149", "DEJ2000": "-50.58425443405" },
        'GD71':      { "_V": "VizieR", "RA_ICRS": "088.11543730306", "DE_ICRS": "+15.88623931523", "Source": "3348071631670500736", "e_RA_ICRS": "0.0445", "e_DE_ICRS": "0.0354", "Plx": "19.5638", "e_Plx": "0.0551", "PM": "189.215", "pmRA": "76.728", "e_pmRA": "0.053", "pmDE": "-172.960", "e_pmDE": "0.038", "RUWE": "0.970", "FG": "118860.92961", "e_FG": "8.6524e+01", "Gmag": "12.999769", "FBP": "9.8669e+04", "e_FBP": "1.2597e+02", "BPmag": "12.853095", "FRP": "3.7757e+04", "e_FRP": "3.4572e+01", "RPmag": "13.305408", "BP-RP": "-0.452312", "RV": " ", "e_RV": " ", "Vbroad": " ", "GRVSmag": " ", "QSO": "0", "Gal": "0", "NSS": "0", "XPcont": "1", "XPsamp": "1", "RVS": "0", "EpochPh": "0", "EpochRV": "0", "MCMCGSP": "0", "MCMCMSC": "0", "And": "0", "Teff": " ", "logg": " ", "[Fe/H]": " ", "Dist": " ", "A0": " ", "HIP": " ", "PS1": "127060881152804422", "SDSS13": " ", "SKYM2": " ", "TYC2": " ", "URAT1": "URAT1-530060563", "AllWISE": "J055227.67+155311.3", "APASS9": "24233759", "GSC23": "N9KN012252", "RAVE5": " ", "2MASS": "05522761+1553137", "RAVE6": " ", "RAJ2000": "88.11508274659", "DEJ2000": "15.88700802429" },
        'EG21':      { "_V": "VizieR", "RA_ICRS": "047.62973152269", "DE_ICRS": "-68.60139793609", "Source": "4646535078125821568", "e_RA_ICRS": "0.0287", "e_DE_ICRS": "0.0254", "Plx": "96.1834", "e_Plx": "0.0289", "PM": "110.596", "pmRA": "39.668", "e_pmRA": "0.036", "pmDE": "-103.237", "e_pmDE": "0.032", "RUWE": "1.127", "FG": "514178.21132", "e_FG": "2.2762e+02", "Gmag": "11.409583", "FBP": "3.9130e+05", "e_FBP": "3.6506e+02", "BPmag": "11.357257", "FRP": "1.9332e+05", "e_FRP": "6.7395e+01", "RPmag": "11.532214", "BP-RP": "-0.174957", "RV": " ", "e_RV": " ", "Vbroad": " ", "GRVSmag": " ", "QSO": "0", "Gal": "0", "NSS": "0", "XPcont": "1", "XPsamp": "1", "RVS": "0", "EpochPh": "0", "EpochRV": "0", "MCMCGSP": "0", "MCMCMSC": "0", "And": "0", "Teff": " ", "logg": " ", "[Fe/H]": " ", "Dist": " ", "A0": " ", "HIP": "14754", "PS1": " ", "SDSS13": " ", "SKYM2": "499272038", "TYC2": "9145-601-1", "URAT1": " ", "AllWISE": "J031031.09-683604.4", "APASS9": "36189953", "GSC23": "S1KQ000128", "RAVE5": " ", "2MASS": "03103100-6836032", "RAVE6": " ", "RAJ2000": "47.62924831674", "DEJ2000": "-68.60093910233" }}

# Coordinates for the affine transformation between GHOST and telescope focal planes
IFU1_CIJ1 = [0.17, 0.979524, -0.000805, -0.15, -0.000161, 0.981456]
IFU2_CIJ1 = [0.23, 0.983066, 0.000644, 0.03, 0.002093, 0.981295]

# Convert from GHOST focal plane to telescope focal plane
def ghost2tel(cij, fpx, fpy):
  a = cij[0]
  b = cij[1]
  c = cij[2]
  d = cij[3]
  e = cij[4]
  f = cij[5]
  g = c/f;

  # Reverse coord transform to get from instrument focalplane to telescope focalplane coord
  x = (fpx - g*fpy + g*d - a)/ (b - g*e)
  y = (fpy - d - e*x)/f

  return (x, y)

# Convert from telescope focal plane to GHOST focal plane
def tel2ghost(cij, x, y):
  # Apply affine transform to correct for alignment and scaling - see astFitij
  fpx = cij[0] + x*cij[1] + y*cij[2]
  fpy = cij[3] + x*cij[4] + y*cij[5]
  return (fpx, fpy)

# Plot the positions of all fibers for a given IFU
# and calculate stats on distances between adjacent fibers
def process_fibers(w, o, header, n, template, color, label, desc, ndist):
    xys = []
    for i in range(n):
        rkey = (template + 'R') % i
        dkey = (template + 'D') % i
        r = f[0].header[rkey]
        d = f[0].header[dkey]
        x, y = np.array(w.wcs_world2pix(r, d, 1)) - o
        if desc[:5] == 'IFU1':
            xys.append(tel2ghost(IFU1_CIJ1, x, y))
        else:
            xys.append(tel2ghost(IFU2_CIJ1, x, y))
        ax.scatter(x, y, marker='x', color=color)
        if ((i == 0) and label):
            ax.annotate(label, [x, y])
    xys = np.array(xys)
    # Calculate the distances between fibers
    # Throw away duplicates and zeros
    d = np.tril(euclidean_distances(xys)).flatten()
    args = d > 0
    # Only keep the number requested - these will be the
    # pairs of closest fibers
    d = np.sort(d[args])[:ndist]
    print(desc, f'n={ndist:d}, mean={d.mean():.4f}, stdev={d.std():.4f} arcsec')

for fn in sys.argv[1:]:
    with fits.open(fn) as f:
        # Get the target name
        target_name = f[0].header['OBJECT']
        print('Target is', target_name)

        # Get the observation time
        tobs_str = f[0].header['DATE'] + 'T' + f[0].header['UTSTART']
        tobs = Time(tobs_str)
        print('Time is', tobs)
        print('--')

        # Get the telescope coordinates from the FITS header
        # These are J2000 and need to have proper motion applied
        dec = float(f[0].header['DEC'])
        telcoord = SkyCoord(ra = float(f[0].header['RA']) * u.deg,
                     dec = dec * u.deg,
                     distance = Distance(parallax = float(f[0].header['PARALLAX']) * u.arcsec),
                     pm_ra_cosdec = float(f[0].header['PMRA']) * u.arcsec/u.yr * 15 * math.cos(math.radians(dec)),
                     pm_dec = float(f[0].header['PMDEC']) * u.arcsec/u.yr,
                     obstime = Time(2000.0, format='jyear'))
        # Apply proper motion
        telcoord_obs = telcoord.apply_space_motion(tobs)
        # Grab just the ra and dec to plot later
        telpos = np.array([telcoord_obs.ra.to_value(), telcoord_obs.dec.to_value()])
        print('Target coordinates (pm applied)', telpos)

        # Do we have GAIA coordinates for this object?
        if target_name in objects:
            target = objects[target_name]
            catcoord = SkyCoord(ra=float(target['RA_ICRS']) * u.deg,
                         dec=float(target['DE_ICRS']) * u.deg,
                         distance=Distance(parallax=float(target['Plx']) * u.mas),
                         pm_ra_cosdec=float(target['pmRA']) * u.mas/u.yr,
                         pm_dec=float(target['pmDE']) * u.mas/u.yr,
                         obstime=Time(2016.0, format='jyear'))
            # Apply proper motion
            catcoord_obs = catcoord.apply_space_motion(tobs)
            # Grab just the ra and dec to plot later
            catpos = np.array([catcoord_obs.ra.to_value(), catcoord_obs.dec.to_value()])
            print('Catalog coordinates (pm applied)', catpos)
        else:
            catpos = None

        # Create a new WCS based on the position of the center of the standard res bundle in IFU1
        w = WCS(naxis = 2)

        # The GHOST focal plane positions of the two IFUs
        ifu1xy = np.array([f[0].header['IFU1X'] - f[0].header['IFU1GDX'], f[0].header['IFU1Y'] - f[0].header['IFU1GDY']])
        ifu2xy = np.array([f[0].header['IFU2X'] - f[0].header['IFU2GDX'], f[0].header['IFU2Y'] - f[0].header['IFU2GDY']])

        # The telescope focal plane positions of the two IFUs
        xy1 = np.array(ghost2tel(IFU1_CIJ1, *ifu1xy))
        xy2 = np.array(ghost2tel(IFU2_CIJ1, *ifu2xy))

        # Use microns in the focal plane as "pixels"
        # And use the position of IFU1 as the fiducial
        w.wcs.crpix = xy1

        # The celestial coordinates of IFU1
        ra = f[0].header['SRIFU1R']
        dec = f[0].header['SRIFU1D']

        # Construct our WCS based on this position
        w.wcs.crval = [ra, dec]

        # degrees per mm (where mm equates to WCS "pixels")
        cdelt = 1.610/3600

        # position angle
        theta = math.radians(f[0].header['PA'] + f[0].header['IAA'])

        w.wcs.cd = [[cdelt * math.cos(theta), -cdelt * math.sin(theta)],
                    [cdelt * math.sin(theta), cdelt * math.cos(theta)]]

        # Standard tangent plane projection
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        w.wcs.cunit = ["deg", "deg"]

        # Conver† the celestial coordinates from the FITS header and from the catalog
        # to xy for plotting
        telxy = np.array(w.wcs_world2pix(*telpos, 1))
        print('Target coordinates in telescope focal plane', telxy)
        if not catpos is None:
            catxy = np.array(w.wcs_world2pix(*catpos, 1))
            print('Catalog coordinates in telescope focal plane', catxy)
            d = np.linalg.norm(telxy - catxy)
            print(f'Distance between target and catalog positions {d:.4e} arcsec')

        '''
        # Some sanity checks
        check = w.wcs_pix2world(0, 0, 1)
        print('check', check)
        check = w.wcs_world2pix(*w.wcs.crval, 1)
        print('check', check)
        '''

        # Create the plot
        fig = plt.figure()
        ax = plt.subplot(projection=w)
        ax.set_xlim([1-175, 1+175])
        ax.set_ylim([1-175, 1+175])
        # Show the patrol area
        circ = mpatches.Circle((1, 1), 274.5/2, edgecolor='blue', facecolor='white')
        ax.add_artist(circ)
        ra_ax = ax.coords[0]
        dec_ax = ax.coords[1]
        ra_ax.set_major_formatter('d.dddddd')
        dec_ax.set_major_formatter('d.dddddd')

        # This is the position of the target based on values in the FITS headers
        ax.scatter(*telxy - (1, 1), marker='.')
        ax.annotate('T', telxy - (1, 1))

        if not catpos is None:
            # This is the position of the target from GAIA DR3
            ax.scatter(*catxy - (1, 1), marker='.')
            ax.annotate('C', catxy - (1, 1))

        # We have FITS headers for the centers of most of the fiber bundles
        ifus = {'SRIFU1': 'IFU1 std res', 'SRIFU2': 'IFU2 std res', 'HRIFU1': 'IFU1 high res', 'SKIFU1': 'IFU1 sky'}
        # Celestial coordinates of each bundle
        coords = {i: [f[0].header[i + 'R'] * u.deg, f[0].header[i + 'D'] * u.deg] for i in ifus}
        # Focal plane coordinates of the center of each bundle
        xys = {i: np.array(w.wcs_world2pix(*coords[i], 1)) for i in ifus}
        # Distance between the target and the center of each bundle
        dists = {i: np.linalg.norm(telxy - xys[i]) for i in ifus}
        # The closest bundle
        arg = min(dists, key=dists.get)
        print(f'Target is closest to {ifus[arg]:s}, {dists[arg]:.4e} arcsec')
        print('--')

        # Print the coordinates of the center of each bundle
        for i, v in ifus.items():
            print(v, coords[i])
            print(v, xys[i])

        # Plot the centers of the bundles
        for i in ifus:
            ax.scatter(*xys[i] - (1, 1), marker='o')

        # Now plot all the fibers
        print("--")
        print("Distances between adjacent fibers in GHOST focal plane")
        print("--")
        process_fibers(w, (1, 1), f[0].header, 7, 'I1SSC%01d', 'r', 'SR', 'IFU1 std res', 12)
        process_fibers(w, (1, 1), f[0].header, 19, 'I1HSC%02d', 'g', 'HR', 'IFU1 high res', 42)
        process_fibers(w, (1, 1), f[0].header, 3, 'I1SSK%01d', 'b', 'SKY', 'IFU1 sky', 3)
        process_fibers(w, (1, 1), f[0].header, 6, 'I1HGD%01d', 'k', None, 'IFU1 high res guide', 6)
        process_fibers(w, (1, 1), f[0].header, 6, 'I1SGD%01d', 'k', None, 'IFU1 std res guide', 6)

        process_fibers(w, (1, 1), f[0].header, 6, 'I2SGD%01d', 'k', None, 'IFU2 std res guide', 6)
        process_fibers(w, (1, 1), f[0].header, 7, 'I2SSC%01d', 'r', 'SR', 'IFU2 std res', 12)
        process_fibers(w, (1, 1), f[0].header, 7, 'I2HSK%01d', 'b', 'SKY', 'IFU2 high res sky', 12)

plt.show()
