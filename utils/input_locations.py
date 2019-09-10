""" 
Python file containing the location of all necessary files to run
the various util scripts that either test aspect of the pipeline, 
fit model parameters or inspect reductions.

This is designed to include a function input to determine which user is
running this. Perhaps a better option will present itself in the future....
"""
import ghostdr.ghost.lookups as lookups
import ghostdr.ghost.lookups.polyfit_dict as polyfit_dict
import astropy.io.fits as pyfits
import os
import numpy as np


class Files():
    """ class that must be ran by each script to define file locations"""

    def __init__(self, user = 'joao', mode='high', cam='red'):
        """ 
        Initialisation function for this class. TThis class has minimal
        error checking and is designed to have all file locations under one
        place, instead of individually in each scrip as it was before. 
        
        Each user should put a different locations depending on their system
        under the user optional part.

        Some things should be common for all users (e.g., ghostdr lookups) and
        are defined as such.
        
        Attributes
        ----------

        user: str, optional
             User that is running the scripts. This should be a command line
             input
        
        mode: str, optional
             Spectrograph mode

        cam: str, optional
             Spectrograph arm

        """

        self.user = user
        self.mode = mode
        self.cam = cam
        # This class has minimal error checking and is designed to have all
        # file locations under one place, instead of individually in each script
        
        # This directory indicates where the ghostdr lookups are.
        # We are using the module path for this, so should be the same for
        # everyone. 
        self.lookups_path = os.path.dirname(os.path.abspath(lookups.__file__))
        self.polyfit_lookups_path = self.lookups_path + '/Polyfit/'
        # Text files containing arc lines, both for all lines and Ar only
        self.arclinefile = self.lookups_path + '/' + lookups.line_list

        # Load all the default parameter files from the lookups. These
        # may need to be overwritten for each individual test with fitted
        # models in the reduced flat or arc.
        self.xmod_location = [value for key, value in
               polyfit_dict.xmod_dict.items()
               if cam in key.lower() and mode in key.lower()][0]
        self.xparams = pyfits.getdata(self.polyfit_lookups_path + self.xmod_location)

        self.wavemod_location = [value for key, value in
               polyfit_dict.wavemod_dict.items()
               if cam in key.lower() and mode in key.lower()][0]
        self.waveparams = pyfits.getdata(self.polyfit_lookups_path + self.wavemod_location)

        self.rotmod_location = [value for key, value in
               polyfit_dict.rotmod_dict.items()
               if cam in key.lower() and mode in key.lower()][0]
        self.rotparams = pyfits.getdata(self.polyfit_lookups_path + self.rotmod_location)

        self.specmod_location = [value for key, value in
               polyfit_dict.specmod_dict.items()
               if cam in key.lower() and mode in key.lower()][0]
        self.specparams = pyfits.getdata(self.polyfit_lookups_path + self.specmod_location)

        self.spatmod_location = [value for key, value in
               polyfit_dict.spatmod_dict.items()
               if cam in key.lower() and mode in key.lower()][0]
        self.spatparams = pyfits.getdata(self.polyfit_lookups_path + self.spatmod_location)

        if self.user=='joao':
            # This is the directory containing the raw (and recently reduced)
            # files.
            self.basedir = '/home/jbento/code/GHOSTDR/simulator/pyghost/output/full_reduction/'
            
            self.arclinefile_ar_only = '/home/jbento/code/GHOSTDR/simulator/pyghost/pyghost/data/mnras_ar_only.txt'

            # Now define locations of actual images for fitting or visualisation.
            # correct)
            self.arc_image_file = self.basedir + "intermediates.Et18qyON6m/arcBefore95_"+self.mode+"_MEF_1x1_"+self.cam+"1_tiled.fits"

            self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
            self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/flat95_'+self.mode+'_1_MEF_1x1_'+self.cam+'1_flat.fits'
            self.flat_image_file = self.basedir + 'calibrations/processed_flat/flat95_'+self.mode+'_1_MEF_1x1_'+self.cam+'1_flat.fits'
            self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
            
            self.slit_flat_image = self.basedir + 'calibrations/processed_slitflat/flat95_'+self.mode+'_2_MEF_2x2_slit_slitflat.fits'
            self.slit_arc_image = self.basedir + 'calibrations/processed_slit/arcBefore95_'+self.mode+'_MEF_2x2_slit_slit.fits'
        elif self.user=='mike':
            # This is the directory containing the raw (and recently reduced)
            # files.
            self.basedir = '/Users/mireland/data/ghost/dhs_testdata_10jul/'
            self.basedir = '/Users/mireland/data/ghost/testdata-clean-190814/'
            
            self.arclinefile_ar_only = '/Users/mireland/python/GHOSTDR/simulator/pyghost/pyghost/data/mnras_ar_only.txt'
            self.default_xmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/high/161120/xmod.fits'
            self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_Xmod.fits'
            self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_Xmod_reversed.fits'
            self.default_wmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/high/161120/wavemod.fits'

            # Now define locations of actual images for fitting or visualisation.
            self.flat_image_file = self.basedir + "flat_processed.fits"
            self.flat_image_file = self.basedir + "processed_flat/flat95_high_1_MEF_1x1_blue1_flat.fits"
            
            
            self.arc_image_file = self.basedir + "intermediates.Et18qyON6m/arcBefore95_"+self.mode+"_MEF_1x1_"+self.cam+"1_tiled.fits"

            self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
            self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/flat95_'+self.mode+'_1_MEF_1x1_'+self.cam+'1_flat.fits'
            #self.flat_image_file = self.basedir + 'calibrations/processed_flat/flat95_'+self.mode+'_1_MEF_1x1_'+self.cam+'1_flat.fits'
            self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
            
            self.slit_flat_image = self.basedir + 'calibrations/processed_slitflat/flat95_'+self.mode+'_2_MEF_2x2_slit_slitflat.fits'
            self.slit_arc_image = self.basedir + 'calibrations/processed_slit/arcBefore95_'+self.mode+'_MEF_2x2_slit_slit.fits'
        else:
            print('Invalid user, try again.')
        
    def thar_spectrum(self, linefile):
        """Calculates a ThAr spectrum. Note that the flux scaling here is roughly
        correct for the lamp with no neutral density. From the simulator.

        Parameters
        ----------

        linefile: string
            Path to the arcline file

        Returns
        -------
        wave, flux: ThAr spectrum (wavelength in um, flux in photons/s?)
        """

        thar = np.loadtxt(linefile, usecols=[0, 1, 2])
        # Create a fixed wavelength scale evenly spaced in log.
        thar_wave = 3600 * np.exp(np.arange(5e5) / 5e5)
        thar_flux = np.zeros(int(5e5))
        # NB This is *not* perfect: we just find the nearest index of the
        # wavelength scale corresponding to each Th/Ar line.
        wave_ix = (np.log(thar[:, 1] / 3600) * 5e5).astype(int)                     
        wave_ix = np.minimum(np.maximum(wave_ix, 0), 5e5 - 1).astype(int)           
        thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))                        
        thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                                mode='same')

        thar_flux /= np.max(thar_flux) / 3e5
        return np.array([thar_wave, thar_flux])
