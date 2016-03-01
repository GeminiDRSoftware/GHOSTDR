class GHOST(DataClassification):
    name="GHOST"
    usage = '''
        Applies to all data from the GHOST instruments in any mode.
        '''
        
    parent = "GEMINI"
    requirement = PHU(INSTRUME='GHOST')

newtypes.append(GHOST())

