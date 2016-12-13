class GHOST_ARC(DataClassification):
    name="GHOST_ARC"
    usage = """
        Applies to all arc datasets from the GHOST instruments
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & PHU(OBSTYPE="ARC")

newtypes.append(GHOST_ARC())
