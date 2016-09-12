class GHOST_DARK(DataClassification):
    name="GHOST_DARK"
    usage = """
        Applies to all flat datasets from the GHOST instruments
        """
    parent = "GHOST"
    requirement = ISCLASS("GHOST") & PHU(OBSTYPE="DARK")

newtypes.append(GHOST_DARK())
