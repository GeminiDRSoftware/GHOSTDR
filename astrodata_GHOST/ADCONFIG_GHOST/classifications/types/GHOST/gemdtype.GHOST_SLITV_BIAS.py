class GHOST_SLITV_BIAS(DataClassification):
    name="GHOST_SLITV_BIAS"
    usage = """
        Applies to all bias datasets from the GHOST slit viewer
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(OBSTYPE="BIAS")

newtypes.append(GHOST_SLITV_BIAS())
