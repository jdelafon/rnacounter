
Usage:
======
See "rnacounter --help" and the tutorial at [...<bbcflib tutorials>].

Installation:
=============

Manually:

    sudo python setup.py install

or with easy_install:

    sudo easy_install rnacounter

or better yet, with pip:

    sudo pip install rnacounter

Dependencies:
=============

* Python 2.7
* setuptools 7.0+

Tests run with the library versions below, but may work with earlier versions.

* numpy 1.6.2+
* docopt 0.6.1+
* cython 0.18+


Content:
========
* rnacounter
The main executable.

* rnacounter.pyx:
Cython version, to be compiled (see setup.py).

* tests/test_rnacounter.py:
Unit tests, run with "nosetests test_rnacounter.py".

* testfiles/:
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after it.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.


Testing:
=========
Unit tests in folder tests/

Example::

    rnacounter testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf

The BAM contains 4041 reads all aligning perfectly on Gapdh (ENSMUSG00000057666) exons,
mostly on ENSMUSE00000487077 but also ENSMUSE00000751942 and ENSMUSE00000886744.
Nothing on other exons, which makes it a good example of badly conditioned input data...

The least squares method returns counts on the following transcripts:
ENSMUST00000117757, ENSMUST00000118875, ENSMUST00000147954
and nothing on ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588 .

Returns a count of 2459.62 (1091.71 RPK) for the gene.

Reports to "benchmark.txt".
