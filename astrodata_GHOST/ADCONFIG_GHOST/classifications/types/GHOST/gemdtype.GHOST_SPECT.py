class GHOST_SPECT(DataClassification):
    name="GHOST_SPECT"
    usage = """
        Applies to all spectra datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & PHU(OBSTYPE="OBJECT")

newtypes.append(GHOST_SPECT())
