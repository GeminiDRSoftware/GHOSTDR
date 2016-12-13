class GHOST_SLITV_FLAT(DataClassification):
    name="GHOST_SLITV_FLAT"
    usage = """
        Applies to all flat datasets from the GHOST slit viewing camera
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(OBSTYPE="FLAT")

newtypes.append(GHOST_SLITV_FLAT())
