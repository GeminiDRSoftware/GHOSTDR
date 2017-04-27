class GHOST_SLITV_ARC(DataClassification):
    name="GHOST_SLITV_ARC"
    usage = """
        Applies to all arc datasets from the GHOST slit viewer
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(OBSTYPE="ARC")

newtypes.append(GHOST_SLITV_ARC())
