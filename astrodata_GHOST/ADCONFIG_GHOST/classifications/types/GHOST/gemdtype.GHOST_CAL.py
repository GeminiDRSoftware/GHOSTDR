class GHOST_CAL(DataClassification):
    name="GHOST_CAL"
    usage = """
        Applies to all calibration datasets from the GHOST instrument
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & OR([
        ISCLASS("GHOST_ARC"),
        ISCLASS("GHOST_BIAS"),
        ISCLASS("GHOST_DARK"),
        ISCLASS("GHOST_FLAT")])

newtypes.append(GHOST_CAL())
