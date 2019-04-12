#!/usr/bin/env python3

from setuptools import setup, find_packages
import sys

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

setup(
    name="plotprofile",
    version="0.1",
    description="Quickly generate reaction path profiles.",
    url="https://github.com/eljost/plotprofile",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "pyyaml",
    ],
    entry_points={
        "console_scripts": [
            "plotprofile = plotprofile.main:run",
        ]
    },
)
