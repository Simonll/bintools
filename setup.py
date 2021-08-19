# coding: utf-8

"""
    Execution Framework

    run

    Contact: evol.simon@gmail.com
"""

from setuptools import find_packages
from setuptools import setup

NAME = "run-utils"
VERSION = "1.0.0"
# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

REQUIRES = ["numpy>=1.18.1", "pandas>=1.0.3"]

setup(
    name=NAME,
    version=VERSION,
    description="simon laurin-lemay run-utils",
    author_email="evol.simon@gmail.com",
    url="https://github.com/Simonll/run/src/",
    install_requires=REQUIRES,
    packages=find_packages(where="src"),
    include_package_data=True,
    package_dir={"": "src"},
)
