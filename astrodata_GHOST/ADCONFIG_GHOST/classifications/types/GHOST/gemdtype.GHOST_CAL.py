class GHOST_CAL(DataClassification):
    name="GHOST_CAL"
    usage = """
        Applies to all calibration datasets from the GHOST instrument
        """
    parent = "GHOST_SPECT"
    requirement = ISCLASS("GHOST_SPECT") & OR([
        ISCLASS("GHOST_ARC"),
        ISCLASS("GHOST_BIAS"),
        ISCLASS("GHOST_DARK"),
        ISCLASS("GHOST_FLAT"),
        ISCLASS("GHOST_SKY"),
    ])

newtypes.append(GHOST_CAL())
