.. _tests:

**********
Test Suite
**********

Available Tests
---------------

The full list of available tests is provided in :ref:`test-suite-API`. In short,
we currently provide:

- Regression tests for descriptors and tag sets;
- Unit tests for primitives, where possible;
- An 'all-up', full reduction test set to check the GHOST pipeline end-to-end.

It is important to note that writing a true 'unit' test for a lot of primitives
isn't possible. Many primitives required a large number of previous steps to
have been completed, proper calibration files to be provided, etc. In such
cases, these primitives are only tested in the 'full reduction' test suite.

Missing or Desirable Tests
--------------------------

The 'full reduction' test suite is currently incomplete, as the details
of the last steps of the reduction process are still being fleshed out.

Running the Tests
-----------------

To invoke the full test suite, simply run::

    pytest

in the top-level directory of the :any:`ghostdr` package.

The 'full reduction' test suite can be run solo by invoking::

    pytest -m fullreduction

The full reduction test suite already possesses the test data that it
requires. The :any:`pytest` fixture which controls the full reduction tests
will automatically place the test files into a temporary directory on your
computer, and then remove all files relating to the tests from that
directory once the tests are complete.

.. warning::
    You'll need 10GB of free disk space to be able to run the full reduction
    test suite. The fixture which prepares the temporary folder for the tests
    will abort the test run if it detects you don't have enough space.

If you only want to run the regression and unit tests, you can invoke
pytest thus::

    pytest -m "not fullreduction"

.. _test-suite-API:

Test Suite API
--------------

Descriptors
^^^^^^^^^^^

.. automodule:: ghost_instruments.test.test_ghost_descriptors
    :members:
    :show-inheritance:

Tags
^^^^

.. automodule:: ghost_instruments.test.test_ghost_tags
    :members:
    :show-inheritance:

Primitive unit tests
^^^^^^^^^^^^^^^^^^^^

.. automodule:: ghostdr.ghost.test
    :members:

Global
""""""

.. automodule:: ghostdr.ghost.test.test_ghost
    :members:

Bundles
"""""""

.. automodule:: ghostdr.ghost.test.test_ghost_bundle
    :members:

Slit frames
"""""""""""

.. automodule:: ghostdr.ghost.test.test_ghost_slit
    :members:

Main spectrograph frames
""""""""""""""""""""""""

.. automodule:: ghostdr.ghost.test.test_ghost_spect
    :members:

Full reduction tests
^^^^^^^^^^^^^^^^^^^^

.. automodule:: ghostdr.ghost.recipes.test
    :members:

Bundles
"""""""

.. automodule:: ghostdr.ghost.recipes.test.000_bundles_test
    :members:

Slit bias
"""""""""

.. automodule:: ghostdr.ghost.recipes.test.010_slitBias_test
    :members:

Slit dark
"""""""""

.. automodule:: ghostdr.ghost.recipes.test.020_slitDark_test
    :members:

Slit flat
"""""""""

.. automodule:: ghostdr.ghost.recipes.test.030_slitFlat_test
    :members:

Slit arc
""""""""

.. automodule:: ghostdr.ghost.recipes.test.040_slitArcStandard_test
    :members:

Standard slit observation
"""""""""""""""""""""""""

.. automodule:: ghostdr.ghost.recipes.test.050_slitObj_test
    :members:

Bias
""""

.. automodule:: ghostdr.ghost.recipes.test.110_overscanCorrect_test
    :members:

.. automodule:: ghostdr.ghost.recipes.test.111_masterBias_test
    :members:

Dark
""""

.. automodule:: ghostdr.ghost.recipes.test.129_masterDark_test
    :members:

Flat
""""

.. automodule:: ghostdr.ghost.recipes.test.139_masterFlat_test
    :members:

Arc
"""

.. automodule:: ghostdr.ghost.recipes.test.149_masterArc_test
    :members:
