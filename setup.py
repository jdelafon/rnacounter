
from distutils.core import setup
from Cython.Build import cythonize
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name = "rnacount",
    ext_modules = cythonize("rnacount.pyx",
                            sources=["rnacounter.cc"],
                            language="c++"),
)
