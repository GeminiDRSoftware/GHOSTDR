class GHOST_BLUE(DataClassification):
    name="GHOST_BLUE"
    usage = '''
        Applies to all data from the blue arm of GHOST instruments in any mode.
        '''
        
    parent = "GHOST_SPECT"
    requirement = ISCLASS('GHOST_SPECT') & PHU(CAMERA='BLUE')

newtypes.append(GHOST_BLUE())

