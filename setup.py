from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np

def readme():
    with open('README.rst') as f:
        return f.read()

extensions = [
    Extension("rnacounter/rnacounter",        # .so
              ["rnacounter/rnacounter.pyx", "rnacounter/rnacounter.c"],  # .c
              include_dirs=[np.get_include()],
             )
]

setup(name='rnacounter',
    version='1.0',
    description='Count reads in genomic intervals',
    long_description=readme(),
    classifiers=[
      'Programming Language :: Python :: 2.7',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='rna-seq rnaseq count reads table sequencing genetics bioinformatics',
    url='https://github.com/delafont/rnacounter',
    author='Julien Delafontaine',
    author_email='julien.delafontaine@epfl.ch',
    license='GPL3',
    packages=['rnacounter'],
    zip_safe=False,
    include_package_data=True,
    test_suite='nose.collector',
    install_requires=['cython','numpy','docopt','nose'],
    scripts=['rnacounter/rnacounter'],
    cmdclass = {'build_ext':build_ext},
    ext_modules = cythonize(extensions),
)


