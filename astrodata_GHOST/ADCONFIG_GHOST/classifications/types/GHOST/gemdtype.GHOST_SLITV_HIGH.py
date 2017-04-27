class GHOST_SLITV_HIGH(DataClassification):
    name="GHOST_SLITV_HIGH"
    usage = """
        Applies to all high res datasets from the GHOST slit viewer
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & PHU(SMPNAME="HI_ONLY")

newtypes.append(GHOST_SLITV_HIGH())
