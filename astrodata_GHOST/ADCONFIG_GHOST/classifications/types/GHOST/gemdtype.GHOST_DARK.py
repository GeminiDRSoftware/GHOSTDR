class GHOST_DARK(DataClassification):
    name="GHOST_DARK"
    usage = """
        Applies to all flat datasets from the GHOST instruments
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & PHU(OBSTYPE="DARK")

newtypes.append(GHOST_DARK())
