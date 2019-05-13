.. recipes:

********************
Recipes and Contexts
********************

Contexts
========

``sq`` - Science Quality
------------------------

``qa`` - Quick Assessment
-------------------------

``test`` - Test Suite
---------------------

.. warning::
    This context is not designed for end-user use, or use outside :any:`pytest`.

It was found during the building of the test suite that a lot of primitives
require full data and calibrators in order to be tested correctly. This meant
that it was best to implement a 'full reduction' test suite.

The best way to do this was to include a ``test`` reduction context, to
allow the RecipeSystem to run 'all-up' in a test environment.

See :ref:`tests` for more details.

Recipes
=======
list of recipes for each context

location of recipes and indexes

Technical Flow Charts
=====================
technical flow charts for each recipes

Issues and Limitations
======================
