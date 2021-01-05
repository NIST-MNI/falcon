#!/usr/bin/env python

import os
import sys

from setuptools import setup, find_packages

os.chdir(os.path.dirname(sys.argv[0]) or ".")

print(find_packages())

setup(
    name = "minc2_simple",
#    packages = ["minc2.simple"],
    version="0.2.30",
    description="MINC2 Simple interface using CFFI",
    long_description=open("README.txt", "rt").read(),
    url="https://github.com/vfonov/minc2_simple",
    author="Vladimir S. FONOV",
    author_email="vladimir.fonov@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
        "License :: OSI Approved :: BSD License",
    ],
    packages=find_packages(),
    install_requires=["cffi>=1.0.0","six"],
    setup_requires=["cffi>=1.0.0","six"],
    cffi_modules=[
        "./minc2_simple/minc2_simple_build.py:ffi",
    ],
    scripts=(['../example/python/xfmavg_scipy.py']),
    test_suite="test"
)
