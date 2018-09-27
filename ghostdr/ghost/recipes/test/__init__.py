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
| Slit bias    | 01n_   |
+--------------+--------+
| Slit dark    | 02n_   |
+--------------+--------+
| Slit flat    | 03n_   |
+--------------+--------+
| Slit arc     | 04n_   |
+--------------+--------+
| Slit         | 05n_   |
+--------------+--------+
| Bias         | 11n_   |
+--------------+--------+
| Dark         | 12n_   |
+--------------+--------+
| Flat         | 13n_   |
+--------------+--------+
| Arc          | 14n_   |
+--------------+--------+
| Standard     | 15n_   |
+--------------+--------+
| Science      | 16n_   |
+--------------+--------+


"""