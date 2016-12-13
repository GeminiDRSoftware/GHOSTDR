class GHOST_SLITV_DARK(DataClassification):
    name="GHOST_SLITV_DARK"
    usage = """
        Applies to all flat datasets from the GHOST slit viewing camera
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(OBSTYPE="DARK")

newtypes.append(GHOST_SLITV_DARK())
