import setuptools
from setuptools import find_packages
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    "astropy",
    "ccdproc",
    # "photutils",
    "multiprocess"
    ]

setuptools.setup(
    name="pyfosc", 
    version="2.0.0-alpha",
    author="Yuming Fu",
    author_email="fuympku@outlook.com",
    description="A pipeline toolbox for long-slit spectroscopy data reduction written in Python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rudolffu/pyfosc",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE V3",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements,
    python_requires='>=3.7',
)
