class GHOST_HIGH(DataClassification):
    name="GHOST_HIGH"
    usage = '''
        Applies to all data from GHOST instruments in high-res mode.
        '''
        
    parent = "GHOST_SPECT"
    requirement = (ISCLASS('GHOST_SPECT') | ISCLASS('GHOST_SLITV')) & PHU(SMPNAME='HI_ONLY')

newtypes.append(GHOST_HIGH())

