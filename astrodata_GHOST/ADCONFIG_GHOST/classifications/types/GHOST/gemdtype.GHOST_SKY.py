class GHOST_SKY(DataClassification):
    name="GHOST_SKY"
    usage = """
        Applies to all sky-frame datasets from the GHOST instruments
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & PHU(OBSTYPE="SKY")

newtypes.append(GHOST_SKY())
