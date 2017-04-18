class GHOST_STD(DataClassification):
    name="GHOST_STD"
    usage = '''
        Applies to all data from GHOST instruments in standard-res mode.
        '''
        
    parent = "GHOST_SPECT"
    requirement = (ISCLASS('GHOST_SPECT') | ISCLASS('GHOST_SLITV')) & PHU(SMPNAME='LO_ONLY')

newtypes.append(GHOST_STD())

