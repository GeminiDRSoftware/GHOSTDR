"""
This is the 'all-up' test set for GHOSTDR.

These tests will effectively do a full data reduction of a simulated test
data set, encompassing all arms and resolution modes. The standard recipes
have been broken apart into smaller recipes, so the output of each can be
tested.

In order the force tests to be executed in the correct order, each module
has a numeric prefix (as pytest will execute the test modules in alphanumeric
order):

+--------------+--------+
| Obs. type    | Prefix |
+--------------+--------+
| Bundles      | 00n\_  |
+--------------+--------+
| Slit bias    | 01n\   |
+--------------+--------+
| Slit dark    | 02n\   |
+--------------+--------+
| Slit flat    | 03n\   |
+--------------+--------+
| Slit arc [*]_| 04n\   |
+--------------+--------+
| Slit [*]_    | 05n\   |
+--------------+--------+
| Bias         | 11n\   |
+--------------+--------+
| Dark         | 12n\   |
+--------------+--------+
| Flat         | 13n\   |
+--------------+--------+
| Arc          | 14n\   |
+--------------+--------+
| Standard     | 15n\   |
+--------------+--------+
| Science      | 16n\   |
+--------------+--------+

.. [*] Note that, because slit arcs and the slits associated standards are
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
