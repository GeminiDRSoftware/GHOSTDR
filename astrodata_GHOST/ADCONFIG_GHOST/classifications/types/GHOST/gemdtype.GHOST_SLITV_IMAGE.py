class GHOST_SLITV_IMAGE(DataClassification):
    name="GHOST_SLITV_IMAGE"
    usage="""
        Applies to a GHOST slit-viewing camera image
        """
    parent = "GHOST_SLITV"
    requirement = ISCLASS("GHOST_SLITV") & NOT(ISCLASS("GHOST_SLITV_CAL"))

newtypes.append(GHOST_SLITV_IMAGE())
