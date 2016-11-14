.. primitives:

**********
Primitives
**********

.. module:: astrodata_GHOST.RECIPES_GHOST.primitives.primitives_GHOST
.. autoclass:: GHOSTPrimitives

    **Primitives**

    .. method:: rejectCosmicRays(self, rc)

        Reject cosmic rays from GHOST data.

        .. warning:: This currently does not successfully flag anything in the
                     DQ plane. I'm looking into this. -MCW

        Parameters
        ----------
        rc : dict
            The ReductionContext dictionary that holds the data stream
            processing information.

        Yields
        ------
        rc : dict
            The same ReductionContext dictionary, with any necessary
            alterations.

    .. method:: stackFrames(self, rc)

        Stack GHOST frames using IRAF.

    .. automethod:: standardizeHeaders

.. **********
.. Primitives
.. **********
..
.. Primitive #1  (alphabetical)
.. ============================
..
.. Test text to make sure this file is being found
..
.. Purpose
.. -------
..
.. Inputs and Outputs
.. ------------------
..
.. Input parameters
.. ----------------
..
.. AstroData Type(s)
.. -----------------
..
.. Inheritance and Primitive Set
.. -----------------------------
..
.. Location
.. --------
..
.. Algorithms
.. ----------
..
.. Issues and Limitations
.. ----------------------
..
..
.. Primitive #2
.. ============
..