# coding: utf-8

"""
    bioinformatics tools

    bintools

    Contact: evol.simon@gmail.com
"""

from setuptools import find_packages
from setuptools import setup

NAME = "bintools"
VERSION = "1.0.0"
# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

REQUIRES = ["numpy>=1.19.2", "pandas>=1.3.4"]

setup(
    name=NAME,
    version=VERSION,
    description="simon laurin-lemay bintools",
    author_email="evol.simon@gmail.com",
    url="https://github.com/Simonll/bintools/src/",
    install_requires=REQUIRES,
    packages=find_packages(where="src"),
    include_package_data=True,
    package_dir={"": "src"},
    # test_suite = "test",
)
