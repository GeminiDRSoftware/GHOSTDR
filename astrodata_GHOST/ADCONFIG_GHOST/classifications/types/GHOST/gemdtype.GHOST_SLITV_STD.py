class GHOST_SLITV_STD(DataClassification):
    name="GHOST_SLITV_STD"
    usage = """
        Applies to all std res datasets from the GHOST slit viewer
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(SMPNAME="LO_ONLY")

newtypes.append(GHOST_SLITV_STD())
