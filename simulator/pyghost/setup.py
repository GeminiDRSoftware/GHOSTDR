# coding: utf-8

""" A Python tool for simulating and extracting GHOST data """

import os
import re
import sys

try:
    from setuptools import setup

except ImportError:
    from distutils.core import setup

major, minor1, minor2, release, serial =  sys.version_info

readfile_kwargs = {"encoding": "utf-8"} if major >= 3 else {}

def readfile(filename):
    with open(filename, **readfile_kwargs) as fp:
        contents = fp.read()
    return contents

version_regex = re.compile("__version__ = \"(.*?)\"")
contents = readfile(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "pyghost",
    "__init__.py"))

version = version_regex.findall(contents)[0]

setup(name="pyghost",
      version=version,
      author="Michael J. Ireland",
      author_email="michael.ireland@anu.edu.au",
      packages=["pyghost"],
      url="http://www.github.com/mikeireland/pyghost/",
      license="MIT",
      description="A Python tool for simulating and extracting GHOST data.",
      long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
      install_requires=[
        "requests",
        "requests_futures"
      ]
     )
