from astrodata import astro_data_tag, TagSet, astro_data_descriptor, returns_list
from gemini_instruments.gemini import AstroDataGemini

from . import lookup
from gemini_instruments.common import build_group_id

class AstroDataGhost(AstroDataGemini):
    """
    Class for adding tags and descriptors to GHOST data.
    """

    __keyword_dict = dict(array_section = 'CCDSEC',
                          array_name = 'AMPNAME',
                          overscan_section = 'BIASSEC',
                          res_mode = 'SMPNAME',
                          exposure_time = 'EXPTIME',
                          )

    @staticmethod
    def _matches_data(source):
        """
        Check if data is from GHOST.

        Parameters
        ----------
        source : astrodata.AstroData
            The source file to check.
        """
        return source[0].header.get('INSTRUME', '').upper() == 'GHOST'

    @astro_data_tag
    def _tag_instrument(self):
        """
        Define the minimal tag set for GHOST data.
        """
        return TagSet(['GHOST'])

    @astro_data_tag
    def _tag_bundle(self):
        """
        Define the 'bundled data' tag set for GHOST data.
        """
        # Gets blocked by tags created by split files
        return TagSet(['BUNDLE'])

    @astro_data_tag
    def _tag_bias(self):
        """
        Define the 'bias data' tag set for GHOST data.
        """
        if self.phu.get('OBSTYPE') == 'BIAS':
            return TagSet(['CAL', 'BIAS'])

    @astro_data_tag
    def _tag_dark(self):
        """
        Define the 'dark data' tag set for GHOST data.
        """
        if self.phu.get('OBSTYPE') == 'DARK':
            return TagSet(['CAL', 'DARK'])

    @astro_data_tag
    def _tag_arc(self):
        """
        Define the 'arc data' tag set for GHOST data.
        """
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['CAL', 'ARC'])

    @astro_data_tag
    def _tag_flat(self):
        """
        Define the 'flat data' tag set for GHOST data.
        """
        if self.phu.get('OBSTYPE') == 'FLAT':
            return TagSet(['CAL', 'FLAT'])

    @astro_data_tag
    def _tag_sky(self):
        """
        Define the 'flat data' tag set for GHOST data.
        """
        if self.phu.get('OBSTYPE') == 'SKY':
            return TagSet(['SKY'])

    @astro_data_tag
    def _tag_res(self):
        """
        Define the tagset for GHOST data of different resolutions.
        """
        if self.phu.get('SMPNAME') == 'HI_ONLY':
            return TagSet(['HIGH'])
        else:
            return TagSet(['STD'])

    @astro_data_tag
    def _tag_slitv(self):
        """
        Define the 'slit data' tag set for GHOST data.
        """
        if self.phu.get('CAMERA', '').lower().startswith('slit'):
            return TagSet(['SLITV', 'IMAGE'], blocks=['SPECT', 'BUNDLE'])

    @astro_data_tag
    def _tag_spect(self):
        """
        Define the 'spectrograph data' tag set for GHOST data.
        """
        # Also returns BLUE or RED if the CAMERA keyword is set thus
        if 'CAMERA' in self.phu:
            return TagSet(({self.phu['CAMERA']} & {'BLUE', 'RED'}) | {'SPECT'},
                          blocks=['BUNDLE'])

    @astro_data_tag
    def _status_processed_ghost_cals(self):
        """
        Define the 'processed data' tag set for GHOST data.
        """
        kwords = set(['PRSLITIM', 'PRSLITBI', 'PRSLITDA', 'PRSLITFL',
                      'PRWAVLFT', 'PRPOLYFT'])
        if set(self.phu) & kwords:
            return TagSet(['PROCESSED'])

    @astro_data_tag
    def _tag_binning_mode(self):
        """
        Define the tagset for GHOST data of different binning modes.
        """
        binnings = self.hdr.get('CCDSUM')
        if isinstance(binnings, list):
            if all([x == binnings[0] for x in binnings]):
                return TagSet([binnings[0].replace(' ', 'x', 1)])
            else:
                return TagSet(['NxN'])
        else:
            return TagSet([binnings.replace(' ', 'x', 1)])

    @astro_data_tag
    def _tag_obsclass(self):
        """
        Define the tagset for 'partnerCal' observations.
        """
        if self.phu.get('OBSCLASS') == 'partnerCal':
            return TagSet(['PARTNER_CAL'])

    @astro_data_descriptor
    def amp_read_area(self):
        """
        Returns a list of amplifier read areas, one per extension, made by
        combining the amplifier name and detector section; or, returns a
        string if called on a single-extension slice.

        Returns
        -------
        list/str
            read_area of each extension
        """
        # Note that tiled arrays won't have an array_name, so we'll fake it
        # FIXME correctly fetch keyword for tileArrays primitive
        if self.phu.get('TILEARRY', None) is not None:
            ampname = [0, ]
        else:
            ampname = self.array_name()
        detsec = self.detector_section(pretty=True)
        # Combine the amp name(s) and detector section(s)
        if self.is_single:
            return "'{}':{}".format(ampname,
                        detsec) if ampname and detsec else None
        else:
            return ["'{}':{}".format(a,d) if a is not None and d is not None else None
                    for a,d in zip(ampname, detsec)]

    @astro_data_descriptor
    def arm(self):
        """
        Returns a string indicating whether these data are from the red
        or blue arm of the spectrograph.

        Returns
        -------
        str/None
            Color of the arm (`'blue'`, `'red'`), or `'slitv'` in case of slit
            viewer data. Returns `None` if arm/slit status can't be determined.
        """
        tags = self.tags
        if 'BLUE' in tags:
            return 'blue'
        elif 'RED' in tags:
            return 'red'
        elif 'SLITV' in tags:
            return 'slitv'
        return None

    @astro_data_descriptor
    def calibration_key(self):
        """
        Returns a suitable calibration key for GHOST, which includes the arm.
        """
        return (self.data_label().replace('_stack', ''), self.arm())

    # FIXME Remove once headers corrected
    @astro_data_descriptor
    def central_wavelength(self, asMicrometers=False, asNanometers=False,
                           asAngstroms=False):
        """
        Dummy to work around current Gemini cal_mgr
        """
        val = self.phu.get(self._keyword_for('central_wavelength'), None)

        if val is None:
            if self.arm() == 'red':
                val = 4000. * 10**-10
            elif self.arm() == 'blue':
                val = 6000. * 10**-10


        if asMicrometers:
            val *= 10**6
        elif asNanometers:
            val *= 10**9
        elif asAngstroms:
            val *= 10**10

        return float(val)

    @astro_data_descriptor
    def detector_name(self):
        """
        Returns the detector (CCD) name.
        """
        if self.phu.get('CCDNAME') is not None:
            return self.phu.get('CCDNAME')
        return self.phu.get('DETTYPE')

    @astro_data_descriptor
    def detector_x_bin(self):
        """
        Returns the detector binning in the x-direction.

        Returns
        -------
        int
            The detector binning
        """
        def _get_xbin(b):
            try:
                return int(b.split()[0])
            except (AttributeError, ValueError):
                return None

        binning = self.hdr.get('CCDSUM')
        if self.is_single:
            return _get_xbin(binning)
        else:
            xbin_list = [_get_xbin(b) for b in binning]
            # Check list is single-valued
            return xbin_list[0] if xbin_list == xbin_list[::-1] else None

    @astro_data_descriptor
    def detector_y_bin(self):
        """
        Returns the detector binning in the y-direction.

        Returns
        -------
        int
            The detector binning
        """
        def _get_ybin(b):
            try:
                return int(b.split()[1])
            except (AttributeError, ValueError, IndexError):
                return None

        binning = self.hdr.get('CCDSUM')
        if self.is_single:
            return _get_ybin(binning)
        else:
            ybin_list = [_get_ybin(b) for b in binning]
            # Check list is single-valued
            return ybin_list[0] if ybin_list == ybin_list[::-1] else None

    # FIXME Remove once headers corrected
    @astro_data_descriptor
    def disperser(self, stripID=False, pretty=False):
        """
        Dummy to work around current Gemini cal_mgr
        """
        return "GHOSTDISP"

    # TODO: GHOST descriptor returns no values if data are unprepared

    @astro_data_descriptor
    def exposure_time(self):
        """
        Returns the exposure time in seconds.

        This function extends the standard exposure_time() descriptor
        by allowing the exposure time to exist in the header of the
        first data extension, as well as the PHU. If the exposure time
        exists in both places, the PHU value takes precedence.

        Returns
        -------
        float
            Exposure time.
        """

        exp_time_default = super(AstroDataGhost, self).exposure_time()

        if exp_time_default is None:
            exposure_time = self[0].hdr.get(
                self._keyword_for('exposure_time'),
                -1)
            if exposure_time == -1:
                return None
            return exposure_time

        return exp_time_default


    # The gain() descriptor is inherited from gemini/adclass, and returns
    # the value of the GAIN keyword (as a list if sent a complete AD object,
    # or as a single value if sent a slice). This is what the GHOST version
    # does so one is not needed here.

    @astro_data_descriptor
    def group_id(self):
        """
        Returns a string representing a group of data that are compatible
        with each other.  This is used when stacking, for example.  Each
        instrument, mode of observation, and data type will have its own rules.

        Returns
        -------
        str
            A group ID for compatible data
        """
        tags = self.tags
        if 'DARK' in tags:
            desc_list = ['exposure_time', 'coadds']
        elif 'BIAS' in tags:
            desc_list = []
        else:  # science exposures (and ARCs)
            desc_list = ['observation_id', 'res_mode']
        desc_list.append('arm')

        # never stack frames of mixed binning modes
        desc_list.append('detector_x_bin')
        desc_list.append('detector_y_bin')

        # MCW: We care about the resolution mode EXCEPT for dark and bias
        if 'DARK' not in tags and 'BIAS' not in tags:
            desc_list.append('res_mode')

        # CJS: Generally need to stop FLATs being stacked with science
        additional_item = 'FLAT' if 'FLAT' in tags else None

        return build_group_id(self, desc_list, prettify=[],
                              additional=additional_item)

    @astro_data_descriptor
    def res_mode(self):
        """
        Get the GHOST resolution mode of this dataset

        Returns
        -------
        str/None
            Resolution of the dataset ('high' | 'std'). Returns `None` if
            resolution mode cannot be determined.
        """
        mode = self.phu.get('SMPNAME')
        try:
            if mode.endswith('HI_ONLY'):
                return 'high'
            elif mode.endswith('LO_ONLY'):
                return 'std'
        except:
            pass
        return None

    # TODO: read_noise(): see comments on gain()

    @astro_data_descriptor
    def read_speed_setting(self):
        """
        GHOST does not require a read speed settings - returns `None`

        Returns
        -------
        `None`
        """
        return None

    @astro_data_descriptor
    def want_before_arc(self):
        """
        This is a special descriptor which is being used as a calibration
        system work-around. Outside of active reduction, this descriptor
        should always return None, as the relevant header keyword should only
        exist very briefly during the fetching of bracketed arc files.

        Returns
        -------
        bool or `None`
        """
        want_before = self.phu.get('ARCBEFOR', None)
        if want_before:
            return True
        elif want_before is None:
            return None
        else:
            return False

