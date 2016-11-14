class GHOST_DARK(DataClassification):
    name="GHOST_DARK"
    usage = """
        Applies to all flat datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & PHU(OBSTYPE="ARC")

newtypes.append(GHOST_DARK())
