class GHOST_BIAS(DataClassification):
    name="GHOST_BIAS"
    usage = """
        Applies to all bias datasets from the GHOST instruments
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & PHU(OBSTYPE="BIAS")

newtypes.append(GHOST_BIAS())
