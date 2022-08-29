#!/usr/bin/env python3

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
import astrodata


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
        # models in the reduced flat or arc. By default, use the most recent
        # parameter file.
        self.xmod_location = [value for key, value in
               polyfit_dict.xmod_dict.items()
               if cam in key.lower() and mode in key.lower()][-1]
        self.xparams = astrodata.open(self.polyfit_lookups_path + self.xmod_location)[0].data

        self.wavemod_location = [value for key, value in
               polyfit_dict.wavemod_dict.items()
               if cam in key.lower() and mode in key.lower()][-1]
        self.waveparams = astrodata.open(self.polyfit_lookups_path + self.wavemod_location)[0].data

        self.rotmod_location = [value for key, value in
               polyfit_dict.rotmod_dict.items()
               if cam in key.lower() and mode in key.lower()][-1]
        self.rotparams = astrodata.open(self.polyfit_lookups_path + self.rotmod_location)[0].data

        self.specmod_location = [value for key, value in
               polyfit_dict.specmod_dict.items()
               if cam in key.lower() and mode in key.lower()][-1]
        self.specparams = astrodata.open(self.polyfit_lookups_path + self.specmod_location)[0].data

        self.spatmod_location = [value for key, value in
               polyfit_dict.spatmod_dict.items()
               if cam in key.lower() and mode in key.lower()][-1]
        self.spatparams = astrodata.open(self.polyfit_lookups_path + self.spatmod_location)[0].data
        
        self.slitvmod_location = [value for key, value in
               polyfit_dict.slitvmod_dict.items()
               if mode in key.lower()][-1]
        self.slitvparams = astrodata.open(self.polyfit_lookups_path + self.slitvmod_location).TABLE[0]
        #import pdb; pdb.set_trace()

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
        elif self.user == 'jon':
            # This is the directory containing the raw (and recently reduced)
            # files.
            self.basedir = '/Volumes/ghost_data/mod/'

            self.arclinefile_ar_only = '/Users/jon/gemini/GHOSTDR/simulator/pyghost/data/mnras_ar_only.txt'
            self.arclinefile = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/ThXe.txt'

            # Now define locations of actual images for fitting or visualisation.
            # correct)
            if self.cam == 'blue':
                if self.mode == 'high':
                    self.arc_image_file = self.basedir + "arc300s_hr100004_arraysTiled.fits"
                    #self.arc_image_file = self.basedir + "blue_hg_hr100012_arraysTiled.fits"
                    #self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
                    self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/cont10s_hr100035_flat.fits'
                    self.flat_image_file = self.basedir + 'calibrations/processed_flat/cont10s_hr100035_flat.fits'
                    #self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
                    self.default_wmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/high/220501/wavemod.fits'
                    self.default_xmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/high/220501/xmod.fits'
                else:
                    self.arc_image_file = self.basedir + "arc300s_lr100008_arraysTiled.fits"
                    #self.arc_image_file = self.basedir + "blue_hg_lr100009_arraysTiled.fits"
                    #self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
                    self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/cont10s_lr100032_flat.fits'
                    self.flat_image_file = self.basedir + 'calibrations/processed_flat/cont10s_lr100032_flat.fits'
                    #self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
                    self.default_wmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/std/220501/wavemod.fits'
                    self.default_xmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/std/220501/xmod.fits'
            else:
                if self.mode == 'high':
                    self.arc_image_file = self.basedir + "arc120s_hr110004_arraysTiled.fits"
                    #self.arc_image_file = self.basedir + "red_hg_hr110015_arraysTiled.fits"
                    #self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
                    self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/cont10s_hr210038_flat.fits'
                    self.flat_image_file = self.basedir + 'calibrations/processed_flat/cont10s_hr210038_flat.fits'
                    #self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
                    self.default_wmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/high/220501/wavemod.fits'
                    self.default_xmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/high/220501/xmod.fits'
                else:
                    self.arc_image_file = self.basedir + "arc120s_lr110011_arraysTiled.fits"
                    #self.arc_image_file = self.basedir + "red_hg_lr110012_arraysTiled.fits"
                    #self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'
                    self.flat_reduced_file = self.basedir + 'calibrations/processed_flat/cont10s_lr210035_flat.fits'
                    self.flat_image_file = self.basedir + 'calibrations/processed_flat/cont10s_lr210035_flat.fits'
                    #self.science_file = self.basedir + 'obj95_0.5_'+self.mode+'_MEF_1x1_'+self.cam+'1_extractedProfile.fits'
                    self.default_wmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/std/220501/wavemod.fits'
                    self.default_xmod = '/Users/jon/gemini/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/std/220501/xmod.fits'
            if self.mode == 'high':
                self.slit_flat_image = self.basedir + 'calibrations/processed_slitflat/sv-20220507.205620-20000_slitflat.fits'
                self.slit_arc_image = self.basedir + 'calibrations/processed_slit/sv-20220507.205620-arc_slit.fits'
            else:
                self.slit_flat_image = self.basedir + 'calibrations/processed_slitflat/sv-20220507.205923-20000_slitflat.fits'
                self.slit_arc_image = self.basedir + 'calibrations/processed_slit/sv-20220507.205923-arc_slit.fits'
        elif self.user=='mike':
            # This is the directory containing the raw (and recently reduced)
            # files.
            self.basedir = '/Users/mireland/data/ghost/dhs_testdata_10jul/'
            self.basedir = '/Users/mireland/data/ghost/testdata-clean-190814/'
            self.basedir = '/Users/mireland/data/ghost/reduced_jun22/calibrations/'
            
            self.arclinefile_ar_only = '/Users/mireland/python/GHOSTDR/simulator/pyghost/pyghost/data/mnras_ar_only.txt'

            self.arc_image_file = self.basedir + "intermediates.Et18qyON6m/arcBefore95_"+self.mode+"_MEF_1x1_"+self.cam+"1_tiled.fits"
            self.arc_reduced_file = self.basedir + 'calibrations/processed_arc/arcBefore95_'+self.mode+'_MEF_1x1_'+self.cam+'1_arc.fits'

            # Now define locations of actual images for fitting or visualisation.
            self.flat_image_file = self.basedir + "flat_processed.fits"
            if self.cam == 'blue':
                self.default_wmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/std/161120/wavemod.fits'
                #self.flat_image_file = self.basedir + "processed_flat/flat95_high_1_MEF_1x1_blue1_flat.fits"
                self.flat_image_file = '/Users/mireland/data/ghost/2October2019/BLUE/calibrations/processed_flat/cont_comb04_flat.fits'
                #191105: Edited line below.
                self.flat_image_file = '/Users/mireland/data/ghost/2019-11-01/BLUE/cont_comb01_flat.fits'
                self.default_xmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/std/161120/xmod.fits'
                self.arc_reduced_file = '/Users/mireland/data/ghost/2October2019/BLUE/calibrations/processed_arc/arc_comb00_arc.fits'
                self.arc_image_file = '/Users/mireland/data/ghost/2October2019/BLUE/arc_comb00_arraysTiled.fits'
                #191105: Edited line below.
                self.arc_image_file = '/Users/mireland/data/ghost/2019-11-04/BLUE/HG/hg_1s00044_arraysTiled.fits'
                self.arclinefile = '/Users/mireland/python/GHOSTDR/utils/Hg.txt'
                self.arc_image_file = '/Users/mireland/data/ghost/2019-11-01/BLUE/arc_comb01_arraysTiled.fits'
                self.arclinefile = '/Users/mireland/python/GHOSTDR/utils/ThXe.txt'
                #191113
                self.arc_image_file = '/Users/mireland/data/ghost/2019-11-08/BLUE/arc_comb01_arraysTiled.fits'
                self.flat_image_file = '/Users/mireland/data/ghost/2019-11-08/BLUE/cont_comb02_flat.fits'
                self.arclinefile = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/Xe.txt'
                #Hg file
                #self.arc_image_file = '/Users/mireland/data/ghost/20sept2019/MONOCHROME/test_Hg_2s_blue_arraysTiled.fits'
                #self.arclinefile = '/Users/mireland/python/GHOSTDR/utils/Hg.txt'
                #On-sky commissioning data!
                self.flat_image_file = self.basedir + 'processed_flat/flats_hr_b6_r6_s01_1x1_20220624_1x1_blue5_flat.fits'
                self.slit_flat_image = self.basedir + 'processed_slitflat/flats_hr_b6_r6_s01_1x1_20220624_2x2_slit_slitflat.fits'
                self.arc_reduced_file = self.basedir + 'processed_arc/arcs_hr_b300_r300_s300_1x1_thxe_2022062_1x1_blue1_arc.fits'
                self.flat_image_file = self.basedir + 'processed_flat/flats_sr_b6_r6_s01_1x1_20220624_1x1_blue5_flat.fits'
                self.slit_flat_image = self.basedir + 'processed_slitflat/flats_sr_b6_r6_s01_1x1_20220624_2x2_slit_slitflat.fits'
                self.arc_reduced_file = self.basedir + 'processed_arc/arcs_sr_b300_r300_s300_1x1_take2_202206_1x1_blue1_arc.fits'
            else:
                self.default_wmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/std/161120/wavemod.fits'
                self.flat_image_file = self.basedir + "processed_flat/flat95_std_1_MEF_1x1_red1_flat.fits"
                self.default_xmod = '/Users/mireland/python/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/std/161120/xmod.fits'
                #191105: Edited below
                self.flat_image_file = "/Users/mireland/data/ghost/2019-11-01/RED/cont_comb02_arraysTiled.fits"
                self.arc_image_file = '/Users/mireland/data/ghost/2019-11-01/RED/arc_comb00_arraysTiled.fits'
                self.arclinefile = '/Users/mireland/python/GHOSTDR/utils/ThXe.txt'
                
                self.flat_image_file = self.basedir + 'processed_flat/flats_sr_b6_r6_s01_1x1_20220624_1x1_red5_flat.fits'
                self.slit_flat_image = self.basedir + 'processed_slitflat/flats_sr_b6_r6_s01_1x1_20220624_2x2_slit_slitflat.fits'
                self.arc_reduced_file = self.basedir + 'processed_arc/arcs_sr_b300_r300_s300_1x1_take2_202206_1x1_red1_arc.fits'
      


        elif self.user=='hayescr':
            # This is the directory containing the raw (and recently reduced)
            # files.
            self.basedir = '/Users/hayesc/research/commissioning/live_reduce/calibrations/'
            
            if (self.cam == 'blue') and (self.mode == 'high'):

                self.arclinefile = '/Users/hayesc/research/GHOSTDR/ghostdr/ghost/lookups/Polyfit/ghost_thar_linelist_20220718.txt '

                self.flat_image_file = self.basedir + 'processed_flat/flat_1x1_hr_br6s02_20220630_1x1_blue1_flat.fits'

                self.slit_flat_image = self.basedir + 'processed_slitflat/flat_1x1_hr_br6s02_20220630_2x2_slit_slitflat.fits'

                self.arc_reduced_file = self.basedir + 'processed_arc/arc_hr_1x1_brs900_morning_20220630_1x1_blue1_arc.fits'

            elif (self.cam == 'blue') and (self.mode == 'std'):

                self.arclinefile = '/Users/hayesc/research/GHOSTDR/ghostdr/ghost/lookups/Polyfit/ghost_thar_linelist_20220718.txt '

                self.flat_image_file = self.basedir + 'processed_flat/flat_1x1_sr_br6s02_set2_20220630_1x1_blue1_flat.fits'

                self.slit_flat_image = self.basedir + 'processed_slitflat/flat_1x1_sr_br6s02_set2_20220630_2x2_slit_slitflat.fits'

                self.arc_reduced_file = self.basedir + 'processed_arc/arc_1x1_sr_br300s150_20220630_1x1_blue1_arc.fits'

            elif (self.cam == 'red') and (self.mode == 'std'):

                self.arclinefile = '/Users/hayesc/research/GHOSTDR/ghostdr/ghost/lookups/Polyfit/ghost_thar_linelist_20220718.txt '

                self.flat_image_file = self.basedir + 'processed_flat/flat_1x1_sr_br6s02_set2_20220630_1x1_red1_flat.fits'

                self.slit_flat_image = self.basedir + 'processed_slitflat/flat_1x1_sr_br6s02_set2_20220630_2x2_slit_slitflat.fits'

                self.arc_reduced_file = self.basedir + 'processed_arc/arc_1x1_sr_br300s150_20220630_1x1_red1_arc.fits'

                self.slit_arc_image = self.basedir + 'processed_slit/arc_1x1_sr_br300s150_20220630_2x2_slit_slit.fits'

                self.arc_image_file = self.basedir + '../arc_1x1_sr_br300s150_20220630_1x1_red1_arraysTiled.fits'


            else:
                self.arclinefile = '/Users/hayesc/research/GHOSTDR/ghostdr/ghost/lookups/Polyfit/ghost_thar_linelist_20220718.txt '

                self.flat_image_file = self.basedir + 'processed_flat/flat_1x1_hr_br6s02_20220630_1x1_red1_flat.fits'

                self.slit_flat_image = self.basedir + 'processed_slitflat/flat_1x1_hr_br6s02_20220630_2x2_slit_slitflat.fits'

                self.arc_reduced_file = self.basedir + 'processed_arc/arc_hr_1x1_brs900_morning_20220630_1x1_red1_arc.fits'

              
            
            #Override Hacks.
            #self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_Xmod.fits'
            #self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_Xmod_reversed.fits'
            #self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_red.fits'
            #self.default_xmod = '/Users/mireland/python/GHOSTDR/utils/new_red_reversed.fits'
           
            
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
        thar_wave = 3600 * np.exp(np.arange(6e5) / 5e5)
        thar_flux = np.zeros(int(6e5))
        # NB This is *not* perfect: we just find the nearest index of the
        # wavelength scale corresponding to each Th/Ar line.
        wave_ix = (np.log(thar[:, 1] / 3600) * 5e5).astype(int)                     
        wave_ix = np.minimum(np.maximum(wave_ix, 0), 6e5 - 1).astype(int)           
        thar_flux[wave_ix] = 10**(np.minimum(thar[:, 2], 4))                        
        thar_flux = np.convolve(thar_flux, [0.2, 0.5, 0.9, 1, 0.9, 0.5, 0.2],
                                mode='same')

        thar_flux /= np.max(thar_flux) / 3e5
        return np.array([thar_wave, thar_flux])
