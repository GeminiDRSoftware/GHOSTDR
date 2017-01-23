class GHOST_SLITV_CAL(DataClassification):
    name="GHOST_SLITV_CAL"
    usage = """
        Applies to all calibration datasets from the GHOST instrument slit-viewing camera
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & OR([
        ISCLASS("GHOST_SLITV_BIAS"),
        ISCLASS("GHOST_SLITV_DARK"),
        ISCLASS("GHOST_SLITV_FLAT")])

newtypes.append(GHOST_SLITV_CAL())
