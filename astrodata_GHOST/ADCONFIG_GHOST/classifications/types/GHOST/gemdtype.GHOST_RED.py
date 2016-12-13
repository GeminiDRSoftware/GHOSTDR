class GHOST_RED(DataClassification):
    name="GHOST_RED"
    usage = '''
        Applies to all data from the red arm of GHOST instruments in any mode.
        '''
        
    parent = "GHOST_SPECT"
    requirement = ISCLASS('GHOST_SPECT') & PHU(CAMERA='RED')

newtypes.append(GHOST_RED())

