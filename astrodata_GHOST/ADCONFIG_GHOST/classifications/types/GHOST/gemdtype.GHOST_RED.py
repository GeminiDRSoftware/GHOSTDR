class GHOST_RED(DataClassification):
    name="GHOST_RED"
    usage = '''
        Applies to all data from the red arm of GHOST instruments in any mode.
        '''
        
    parent = "GHOST"
    requirement = ISCLASS('GHOST') & PHU(CAMERA='RED')

newtypes.append(GHOST_RED())

