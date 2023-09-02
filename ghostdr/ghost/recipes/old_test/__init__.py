"""
This is the 'all-up' old_test set for GHOSTDR.

These tests will effectively do a full data reduction of a simulated old_test
data set, encompassing all arms and resolution modes. The standard recipes
have been broken apart into smaller recipes, so the output of each can be
tested.

In order the force tests to be executed in the correct order, each module
has a numeric prefix (as pytest will execute the old_test modules in alphanumeric
order):

+--------------+--------+
| Obs. type    | Prefix |
+--------------+--------+
| Bundles      | 00n\_  |
+--------------+--------+
| Slit bias    | 01n\_  |
+--------------+--------+
| Slit dark    | 02n\_  |
+--------------+--------+
| Slit flat    | 03n\_  |
+--------------+--------+
| Slit arc [1]_| 04n\_  |
+--------------+--------+
| Slit [1]_    | 05n\_  |
+--------------+--------+
| Bias         | 11n\_  |
+--------------+--------+
| Dark         | 12n\_  |
+--------------+--------+
| Flat         | 13n\_  |
+--------------+--------+
| Arc          | 14n\_  |
+--------------+--------+
| Standard     | 15n\_  |
+--------------+--------+
| Science      | 16n\_  |
+--------------+--------+

.. [1] Note that, because slit arcs and the slits associated standards are
       treated identically, they are currently tested together
       under the slit arc heading.

All tests are marked in pytest as 'fullreduction' - to run only the full
reduction tests, invoke pytest as follows::

    pytest -m fullreduction

Alternatively, to skip the full reduction tests, invoke pytest thus::

    pytest -m "not fullreduction"

"""

# Make sure the GHOST instruments package gets registered
import ghost_instruments
import ghostdr
import sqlite3


def get_caldb_contents(dbpath):
    print(dbpath)
    conn = sqlite3.connect(dbpath)
    c = conn.cursor()
    c.execute('SELECT * FROM diskfile')
    return c.fetchall()
