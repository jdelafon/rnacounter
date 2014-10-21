from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np

def readme():
    with open('README.rst') as f:
        return f.read()

extensions = [
    Extension("rnacounter/rnacounter",        # .so
              ["rnacounter/rnacounter.pyx"],  # .c
              include_dirs=[np.get_include()],
             )
]

setup(name='rnacounter',
    version='1.0',
    description='Estimate abundances of genomic features from read densities',
    long_description=readme(),
    classifiers=[
      'Programming Language :: Python',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Natural Language :: English',
      'Operating System :: OS Independent',
      'Development Status :: 5 - Production/Stable',
    ],
    keywords='rna-seq rnaseq count reads table sequencing genetics bioinformatics',
    url='https://github.com/delafont/rnacounter',
    author='Julien Delafontaine',
    author_email='julien.delafontaine@epfl.ch',
    license='GPL-2',
    packages=['rnacounter'],
    zip_safe=False,
    include_package_data=True,
    test_suite='nose.collector',
    install_requires=['cython','numpy','scipy','docopt','nose','pysam'],
    scripts=['rnacounter/rnacounter', 'rnacounter/rnacounter3'],
    cmdclass = {'build_ext':build_ext},
    ext_modules = cythonize(extensions),
)


