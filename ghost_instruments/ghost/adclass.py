from astrodata import astro_data_tag, TagSet, astro_data_descriptor, returns_list
from gemini_instruments.gemini import AstroDataGemini

from . import lookup
from gemini_instruments.common import build_group_id

class AstroDataGhost(AstroDataGemini):

    __keyword_dict = dict(array_section = 'CCDSEC',
                          array_name = 'AMPNAME',
                          overscan_section = 'BIASSEC',
                          res_mode = 'SMPNAME',
                          )

    @staticmethod
    def _matches_data(data_provider):
        return data_provider.phu.get('INSTRUME', '').upper() == 'GHOST'

    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['GHOST'])

    @astro_data_tag
    def _tag_bundle(self):
        # Gets blocked by tags created by split files
        return TagSet(['BUNDLE'])

    @astro_data_tag
    def _tag_bias(self):
        if self.phu.get('OBSTYPE') == 'BIAS':
            return TagSet(['CAL', 'BIAS'])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('OBSTYPE') == 'DARK':
            return TagSet(['CAL', 'DARK'])

    @astro_data_tag
    def _tag_arc(self):
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['CAL', 'ARC'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('OBSTYPE') == 'FLAT':
            return TagSet(['CAL', 'FLAT'])

    @astro_data_tag
    def _tag_sky(self):
        if self.phu.get('OBSTYPE') == 'SKY':
            return TagSet(['SKY'])

    @astro_data_tag
    def _tag_res(self):
        if self.phu.get('SMPNAME') == 'HI_ONLY':
            return TagSet(['HIGH'])
        else:
            return TagSet(['STD'])

    @astro_data_tag
    def _tag_slitv(self):
        if self.phu.get('CCDNAME', '').startswith('Sony-ICX674'):
            return TagSet(['SLITV', 'IMAGE'], blocks=['SPECT', 'BUNDLE'])

    @astro_data_tag
    def _tag_spect(self):
        # Also returns BLUE or RED if the CAMERA keyword is set thus
        if 'CAMERA' in self.phu:
            return TagSet(({self.phu['CAMERA']} & {'BLUE', 'RED'}) | {'SPECT'}, blocks=['BUNDLE'])

    @astro_data_tag
    def _status_processed_ghost_cals(self):
        kwords = set(['PRSLITIM', 'PRSLITBI', 'PRSLITDA', 'PRSLITFL',
                      'PRWAVLFT', 'PRPOLYFT'])
        if set(self.phu.keywords) & kwords:
            return TagSet(['PROCESSED'])

    @astro_data_descriptor
    def amp_read_area(self):
        """
        Returns a list of amplifier read areas, one per extension, made by
        combining the amplifier name and detector section. Or returns a
        string if called on a single-extension slice.

        Returns
        -------
        list/str
            read_area of each extension
        """
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
            Color of the arm ('blue' | 'red')
        """
        tags = self.tags
        if 'BLUE' in tags:
            return 'blue'
        elif 'RED' in tags:
            return 'red'
        elif 'SLITV' in tags:
            return 'slit'
        return None

    @astro_data_descriptor
    def calibration_key(self):
        """
        Returns a suitable calibration key for GHOST, which includes the arm.
        """
        return (self.data_label().replace('_stack', ''), self.arm())

    @astro_data_descriptor
    def detector_x_bin(self):
        """
        Returns the detector binning in the x-direction

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
        Returns the detector binning in the y-direction

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

    # TODO: GHOST descriptor returns no values if data are unprepared
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
            Resolution of the dataset ('high' | 'std')
        """
        mode = self.phu.get('SMPNAME')
        if mode == 'HI_ONLY':
            return 'high'
        elif mode == 'LO_ONLY':
            return 'std'
        return None

    # TODO: read_noise(): see comments on gain()

    @astro_data_descriptor
    def read_speed_setting(self):
        return None

