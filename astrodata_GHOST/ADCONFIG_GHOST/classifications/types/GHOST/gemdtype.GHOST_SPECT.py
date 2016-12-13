class GHOST_SPECT(DataClassification):
    name="GHOST_SPECT"
    usage = """
        Applies to all spectra datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & PHU(DETTYPE="E2V-CCD-231-84")

newtypes.append(GHOST_SPECT())
