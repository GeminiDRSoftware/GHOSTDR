class GHOST_SLITV(DataClassification):
    name="GHOST_SLITV"
    usage = """
        Applies to all slit viewer datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & PHU(CCDNAME="^Sony-ICX674.*$")

newtypes.append(GHOST_SLITV())
