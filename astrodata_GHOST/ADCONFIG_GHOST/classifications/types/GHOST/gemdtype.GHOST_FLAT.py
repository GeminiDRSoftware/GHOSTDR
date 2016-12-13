class GHOST_FLAT(DataClassification):
    name="GHOST_FLAT"
    usage = """
        Applies to all flat datasets from the GHOST instruments
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & PHU(OBSTYPE="FLAT")

newtypes.append(GHOST_FLAT())
