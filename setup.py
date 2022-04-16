from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

ext_modules = [
    Pybind11Extension("_psfmodels", ["src/_psfmodels/pythonBindings.cpp"], cxx_std=17)
]

setup(ext_modules=ext_modules)
