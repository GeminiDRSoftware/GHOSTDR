# This parameter file contains the parameters related to the primitives located
# in the primitives_ghost_spect.py file, in alphabetical order.

from gempy.library import config

from geminidr.core import (
    parameters_ccd, parameters_visualize, parameters_preprocess,
    parameters_spect, parameters_stack)

from astrodata import AstroData as ad

def arcs_valueCheck(value):
    """Validate applyFlatBPMConfig.arcs"""
    return len(value) == 2 and isinstance(
        value[0], str) and isinstance(
        value[1], str)


class addWavelengthSolutionConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_wavelengthAdded",
                          optional=True)
    arc_before = config.ListField("Before arc to use for wavelength solution",
                            (str, ad), None, optional=True, single=True)
    arc_after = config.ListField("After arc to use for wavelength solution",
                            (str, ad), None, optional=True, single=True)


class applyFlatBPMConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_flatBPMApplied",
                          optional=True)
    flat = config.ListField("Flat field to use", (str, ad), None,
                            optional=True, single=True)
    flat_stream = config.Field("Stream to obtain flat field from", str, None,
                               optional=True)
    write_result = config.Field("Write primitive output to disk?", bool, True,
                                optional=True)


class barycentricCorrectConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_barycentricCorrected",
                          optional=True)
    velocity = config.Field("Radial velocity correction", float,
                            None, optional=True)


class clipSigmaBPMConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_sigmaBPMClipped",
                          optional=True)
    sigma = config.Field("Sigma value for clipping", float, 3.0)
    iters = config.Field("Number of clipping iterations", int, None,
                         optional=True)
    bpm_value = config.Field("BPM value to give to clipped pixels", int, 1)


class darkCorrectConfig(parameters_preprocess.darkCorrectConfig):
    def setDefaults(self):
        self.suffix = "_darkCorrected"
        self.dark = None
        self.do_cal = "skip"


class extractProfileConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_extracted",
                          optional=True)
    slit = config.ListField("Slit viewer exposure", (str, ad), None,
                            optional=True, single=True)
    slitflat = config.ListField("Slit viewer flat field", (str, ad), None,
                                optional=True, single=True)
    flat = config.ListField("Flat field", (str, ad), None,
                            optional=True, single=True)
    ifu1 = config.ChoiceField("Status of IFU1", str,
                              allowed={"object": "Pointed at an object",
                                       "sky": "Pointed at sky",
                                       "stowed": "Stowed"},
                              default=None, optional=True)
    ifu2 = config.ChoiceField("Status of IFU1", str,
                              allowed={"object": "Pointed at an object",
                                       "sky": "Pointed at sky",
                                       "stowed": "Stowed"},
                              default=None, optional=True)
    sky_subtract = config.Field("Sky-subtract object spectra?", bool, True)
    flat_correct = config.Field("Flatfield the data?", bool, True)
    snoise = config.RangeField("Fraction of signal to be added to noise estimate for CR flagging",
                               float, 0.1, min=0, max=1)
    sigma = config.RangeField("Number of standard deviations at which to flag pixels",
                              float, 6, min=3)
    weighting = config.ChoiceField("Pixel weighting scheme for extraction", str,
                                   allowed={"uniform": "uniform weighting",
                                            "optimal": "optimal extraction"},
                                   default="optimal")
    tolerance = config.RangeField("Fractional tolerance for convergence",
                                  float, 0.001, min=1e-8, max=0.05)
    apply_centroids = config.Field("Apply slit center-of-light offsets?", bool, False)
    seeing = config.RangeField("FWHM of seeing disc if no processed_slit is "
                               "available", float, None, min=0.2, optional=True)
    write_result = config.Field("Write primitive output to disk?", bool, False)
    debug_cr_map = config.Field("Add CR map to output?", bool, False)
    debug_order = config.RangeField("Order for CR debugging plot", int, None,
                                       min=33, max=97, optional=True)
    debug_pixel = config.RangeField("Pixel for CR debugging plot", int, None,
                                       min=0, max=6144, optional=True)
    debug_timing = config.Field("Output time per order?", bool, False)


class interpolateAndCombineConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_interpdAndCombined",
                          optional=True)
    scale = config.ChoiceField("Output wavelength scale", str, {
        'linear': 'Linear wavelength scale',
        'loglinear': 'Log-linear wavelength scale'
    }, default='loglinear')
    oversample = config.Field("(Approx.) oversampling of output wavelength "
                              "scale", float, 1.0)


class findAperturesConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_aperturesFound",
                          optional=True)
    slitflat = config.Field("Slit viewer flat field",
                            (str, ad),
                            None, optional=True)
    flat = config.ListField("Flat field", (str, ad), None,
                            optional=True, single=True)
    make_pixel_model = config.Field('Add a pixel model to the flat field?',
                                    bool, False)


class fitWavelengthConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_wavelengthFitted",
                          optional=True)
    flat = config.ListField("Flat field", (str, ad), None,
                            optional=True, single=True)
    min_snr = config.RangeField("Minimum S/N for peak detection",
                                float, 20, min=10)
    sigma = config.RangeField("Number of standard deviations for rejecting lines",
                              float, 3, min=1)
    max_iters = config.RangeField("Maximum number of iterations", int, 1, min=1)
    radius = config.RangeField("Matching distance for lines", int, 12, min=2)
    plot1d = config.Field("Produce 1D plots of each order to inspect fit?",
                          bool, False)
    plotrms = config.Field("Produce rms scattergram to inspect fit?",
                           bool, False)
    debug_plot2d = config.Field("Produce 2D plot to inspect fit?", bool, False)


class flatCorrectConfig(config.Config):
    skip = config.Field("No-op this primitive?", bool, False, optional=True)
    suffix = config.Field("Filename suffix", str, "_flatCorrected",
                          optional=True)
    slit = config.Field("Slit viewer exposure", (str, ad), None, optional=True)
    slitflat = config.ListField("Slit viewer flat field", (str, ad), None,
                                optional=True, single=True)
    flat = config.ListField("Processed flat field exposure", (str, ad), None,
                            optional=True, single=True)
    write_result = config.Field("Write primitive output to disk?", bool, True,
                                optional=True)


class formatOutputConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_formattedOutput",
                          optional=True)
    detail = config.ChoiceField("Level of detail", str, {
        'default': 'Default output',
        'processed_image': 'Include processed CCD image',
        'flat_profile': 'Include flat profile',
        'sensitivity_curve': 'Include computed sensitivity curve',
    }, default='default')


class measureBlazeConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_blazeMeasured",
                          optional=True)
    slitflat = config.ListField("Slit viewer flat field", (str, ad), None,
                                optional=True, single=True)


class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    def setDefaults(self):
        self.suffix = '_overscanCorrect'
        self.niterate = 2
        self.high_reject = 3.0
        self.low_reject = 3.0
        self.function = 'poly'
        self.nbiascontam = 4
        self.order = 0


class removeScatteredLightConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_scatteredLightRemoved",
                          optional=True)
    skip = config.Field("Skip removal of scattered light?", bool, True)
    debug_spline_smoothness = config.RangeField(
        "Scaling factor for spline smoothness", float, default=1, min=0.5)
    debug_save_model = config.Field("Attach scattered light model to output?",
                                    bool, False)


class responseCorrectConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_responseCorrected",
                          optional=True)
    standard = config.Field("Standard star (observed)", (str, ad), None,
                       optional=True)
    specphot_file = config.Field("Filename containing spectrophotometry", str, None,
                            optional=True)
    units = config.Field("Units for output spectrum", str, "W m-2 nm-1",
                         check=parameters_spect.flux_units_check)
    write_result = config.Field("Write primitive output to disk?", bool, False)
    order = config.RangeField("Order of polynomial fit to each echelle order", int,
                              1, min=1, max=5)
    debug_plots = config.Field("Show response-fitting plots for each order?",
                               bool, False)


class stackArcsConfig(parameters_stack.core_stacking_config):
    skip = config.Field("No-op this primitive?", bool, False)
    time_delta = config.RangeField("Max. time separating bracketed arcs (seconds)",
                              float, 1200, min=0, optional=True)
    write_result = config.Field("Write primitive output to disk?", bool, True,
                                optional=True)

    def setDefaults(self):
        self.operation = "lmedian"


class standardizeSpectralFormatConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_dragons",
                          optional=True)


write1DSpectraConfig = parameters_spect.write1DSpectraConfig


class createFITSWCSConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_wfits",
                          optional=True)
    iraf = config.Field("Use IRAF (non-standard) format?", bool, False)
    angstroms = config.Field("Write wavelength as Angstroms?", bool, False)


class tileArraysConfig(parameters_visualize.tileArraysConfig):
    def setDefaults(self):
        self.suffix = "_arraysTiled"
