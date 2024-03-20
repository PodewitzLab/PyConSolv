#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import path
from setuptools import setup, find_packages
import sys

min_version = (3, 10)
if sys.version_info < min_version:
    error = """
PyConSolv does not support Python {0}.{1}.
Python {2}.{3} and above is required. Check your Python version like so:
python3 --version
This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:
pip install --upgrade pip
""".format(*(sys.version_info[:2] + min_version))
    sys.exit(error)

here = path.dirname(path.realpath(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'requirements.txt')) as requirements_file:
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]
setup(
    name="PyConSolv",
    version="1.0.6.3",
    description="A package for conformer generation of transition-metal-containing complexes",
    long_description=readme,
    author="Radu Alexandru Talmazan",
    author_email="radu.talmazan@tuwien.ac.at",
    python_requires='>={}'.format('.'.join(str(n) for n in min_version)),
    install_requires=requirements,
    zip_safe=False,
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    packages=['PyConSolv',
              'PyConSolv.interfaces',
              'PyConSolv.misc',
              'PyConSolv.utils'],
    package_data={
        'PyConSolv': [
            'Manual/*'
        ]
    },
    py_modules=[
        "PyConSolv.interfaces",
        "PyConSolv.misc",
        "PyConSolv.utils",
    ],
    entry_points={
        'console_scripts': [
            'pyconsolv = PyConSolv.pyconsolv:main',
        ],
    },
)