
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "rnacount",
    ext_modules = cythonize("rnacount.pyx"),
)
