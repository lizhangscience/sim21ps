#!/usr/bin/evn python3
#
# References:
# [1] Python Pakcaging User Guide: packaging and distributing projects
#     https://packaging.python.org/distributing/
# [2] fg21sim codes
#     https://github.com/liweitianux/fg21sim/setup.py

import os

from setuptools import setup, find_packages

import fg21sim as pkg


def read(fname):
    """
    read files  (to be rewrited)
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name=pkg.__pkgname__,
    version=pkg.__version__,
    description=pkg.__description__,
    long_description=read("README.rst"),
    author=pkg.__author__,
    author_email=pkg.__author_email__,
    url=pkg.__url__,
    license=pkg.__license__,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    packages=find_packages(exclude=["docs", "tests"]),
    scripts=[
        "bin/simAGN",
    ],
    install_requires=[
        "numpy",
        "scipy",
        "matploglib",
        "astropy",
        "healpy",
        "configobj",
    ],
    tests_require=["pytest"],
)
