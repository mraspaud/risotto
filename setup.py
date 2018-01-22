#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Setup stuff
"""

from setuptools import setup
from distutils.extension import Extension
import os

try:
    from Cython.Build import cythonize
    ext = '.pyx'
    USE_CYTHON = True
except ImportError:
    ext = '.c'
    USE_CYTHON = False

extensions = [Extension("risotto._risotto", [os.path.join("risotto/_risotto"+ext)])]

if USE_CYTHON:
    extensions = cythonize(extensions)

setup(
    name='risotto',
    version='0.0.3',
    description="Rice aka USES aka SZIP decoder in Python/Cython",
    author='Martin Raspaud',
    author_email='martin.raspaud@smhi.se',
    url="https://github.com/mraspaud/risotto",
    ext_modules=extensions,
    packages=['risotto'],
    install_requires=['numpy'],
    classifiers=["Development Status :: 4 - Beta",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 " +
                 "or later (GPLv3+)",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python",
                 "Topic :: Scientific/Engineering",
                 "Topic :: System :: Archiving :: Compression"],
)
