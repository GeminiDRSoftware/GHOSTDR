class GHOST_OBJECT(DataClassification):
    name="GHOST_OBJECT"
    usage="""
        Applies to a GHOST science observation
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & NOT(ISCLASS("GHOST_CAL"))

newtypes.append(GHOST_OBJECT())
