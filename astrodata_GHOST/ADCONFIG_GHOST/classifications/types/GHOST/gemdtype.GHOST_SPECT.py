class GHOST_SPECT(DataClassification):
    name="GHOST_SPECT"
    usage = """
        Applies to all spectra datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & OR([
        PHU(OBSTYPE="FLAT"),
        PHU(OBSTYPE="ARC"),
        PHU(OBSTYPE="OBJECT"),
        PHU(OBSTYPE="SKY")])

newtypes.append(GHOST_SPECT())
