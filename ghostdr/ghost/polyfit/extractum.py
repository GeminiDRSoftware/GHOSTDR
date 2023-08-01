import numpy as np
from scipy.ndimage import convolve1d
from scipy import optimize
from amtplotlib import pyplot as plt

from gempy.gemini.library import astrotools as at


class Extractum:
    """
    This is a class that contains information about the data at a single
    wavelength in the echellogram.
    """
    spat_conv_weights = np.array([0.25, .5, .25])

    def __init__(self, phi, data, inv_var=None, mask=None,
                 noise_model=None, pixel=None):
        self.phi = phi
        self.nprof, self.npix = phi.shape
        self.data = data
        self.badpix = np.zeros((self.npix,), dtype=bool) if mask is None else mask
        self.noise_model = noise_model
        self.inv_var = at.divide0(1., noise_model(data)) if inv_var is None and data is not None else inv_var
        self.cr = np.zeros_like(self.badpix)
        self.pixel = pixel

    def fit(self, good=None, coeffs=None):
        if good is None:
            good = ~(self.badpix | self.cr)
        if good.sum() < self.nprof:
            good = None

        # Iterate?

        # Don't use the weights initially because these are smoothed
        # and small CRs get missed
        if coeffs is None:
            result = optimize.lsq_linear(self.phi.T[good], self.data[good],
                                         bounds=(-np.inf, np.inf))
            coeffs = result.x
        model = np.dot(coeffs, self.phi)
        #inv_var_use = convolve1d(at.divide0(1., self.noise_model(model)),
        #                         self.spat_conv_weights)
        inv_var_use = at.divide0(1., self.noise_model(model))
        A = self.phi * np.sqrt(inv_var_use)
        b = self.data * np.sqrt(inv_var_use)
        result = optimize.lsq_linear(A.T[good], b[good], bounds=(-np.inf, np.inf))
        return result.x

    def find_cosmic_rays(self, snoise=0.1, sigma=6, debug=False):
        sigmasq = sigma * sigma

        # Don't try to fit points which aren't expected to have any signal
        good = np.logical_and.reduce([self.inv_var > 0,
                                      self.phi.sum(axis=0) > 0,
                                      ~self.badpix, ~self.cr])

        # If we don't have enough pixels to fit we can't reject
        if good.sum() < self.nprof:
            return []

        try:
            coeffs, _ = optimize.nnls(self.phi.T[good], self.col_data[good])
        except:
            print(f"INITIAL FITTING FAILURE AT PIXEL {self.pixel}")
            print(self.data)
            for p in self.phi:
                print(p)
            print(good)
            raise

        orig_good = good.copy()
        if debug:
            print("\n")
            print("OBJECT PROFILES (phi)")
            print(self.phi)

        model = np.dot(coeffs, self.phi)
        inv_var_use = 1. / np.maximum(self.noise_model(model) + (snoise * model) ** 2, 0)
        nit = 0
        while good.sum() >= self.nprof:
            # CJS problem with using weights and rejecting low pixels
            # Further investigation needed... can we be sure there are
            # no low pixels?
            A = self.phi * np.sqrt(inv_var_use)
            b = self.col_data * np.sqrt(inv_var_use)
            try:
                x, _ = optimize.nnls(A.T[good], b[good])
            except:
                print(f"FITTING FAILURE AT PIXEL {self.pixel} ITERATION {nit}")
                model = np.dot(coeffs, self.phi)
                inv_var_use = convolve1d(at.divide0(1., self.noise_model(model)),
                                         self.spat_conv_weights)
                print(self.data[good])
                for p in self.phi:
                    print(p[good])
                print("-" * 60)
                print(inv_var_use)
                print(model)
                print("-" * 60)
                A = self.phi * np.sqrt(inv_var_use)
                b = self.data * np.sqrt(inv_var_use)
                print(b[good])
                for p in A:
                    print(p[good])
                raise

            model = np.dot(coeffs, self.phi)
            if debug:
                print("RAW COLUMN DATA")
                print(self.data)
                print("RESULT", coeffs)
                print(model)
                print(good)
                print("-" * 60)

            cr_var = convolve1d(np.maximum(self.noise_model(model) + (snoise * model) ** 2, 0),
                                self.spat_conv_weights)
            inv_cr_var = np.where(cr_var > 0, 1. / cr_var, 0)
            deviation = (self.data - model) ** 2 * inv_cr_var  # no sqrt()
            if debug:
                print(np.sqrt(cr_var))
                print(deviation)
                print("-" * 60)
            worst_offender = np.argmax(abs(deviation)[good])
            if abs(deviation[good][worst_offender]) > sigmasq:
                good[np.where(good)[0][worst_offender]] = False
            else:
                break
            nit += 1

        new_bad = np.where(np.logical_and.reduce(
            [self.inv_var > 0, abs(deviation) > sigmasq, orig_good]))[0]
        self.cr[new_bad] = True

        # CJS 20230120: this allows a little extra leeway for vertical CCD bleed
        # limit = ndimage.maximum_filter(y_hat + nsigma * np.sqrt(var_use),
        #                               size=3, mode='constant')
        # limit = y_hat + nsigma * np.sqrt(var_use)
        # new_bad = np.where(col_data > limit)[0]
        # if debug and len(new_bad) > 0:
        if debug:
            print("BAD", new_bad, self.phi.shape)
            # Sky is always last profile in phi?
            summed_fluxes = [np.sum(self.col_data[p > 0]) for p in self.phi[:-1]]
            print("SUMMED FLUXES IN OBJECT FIBRES (not sky subtracted)",
                  summed_fluxes)
            plt.ioff()
            plt.plot(self.col_data, 'k-', label='data')
            # max allowed value
            plt.plot(model + sigma * np.sqrt(cr_var), 'r-', label='limit')
            plt.plot(model, 'b-', label='model')
            for ppp, amp in zip(self.phi, coeffs):
                plt.plot(ppp * amp, ':')
            plt.show()
            plt.ion()

        return new_bad
