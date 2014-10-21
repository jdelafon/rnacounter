
`rnacounter` estimates abundances of genes and their different transcripts
from read densities. Exons and introns can also be quantified.
It requires a genome-level BAM file and a
GTF/GFF file describing the exon structure, such as those provided by Ensembl or GenRep.
The method used is described in [<ref>].

The code is under GPL-2 license.

Usage:
======
See "rnacounter --help" and the tutorial at
`http://bbcf.epfl.ch/bbcflib/tutorial_rnacounter.html`_,
also available in the doc/ folder.

Minimal example::

    rnacounter test.bam test.gtf        # Python 2.7 version
    rnacounter3 test.bam test.gtf       # Python 3 version

Installation:
=============
Manually:

    sudo python setup.py install

or better, with easy_install:

    sudo easy_install rnacounter

or better yet, with pip:

    sudo pip install rnacounter

It installs as a standard Python library but includes the executable
and puts it somewhere in your $PATH. Dependencies should be added
automatically.

Use "python3", "easy_install3", "pip3" respectively to run it with the latest version of Python.
The code is fully compatible with Python 2.7, but uses Python3's syntax.
In particular, things that became iterators in Python3 will be treated as lists
if Python2.7 is used (which may increase speed at the expense of the memory used).

Dependencies:
=============
Tests run with the library versions below, but may work with earlier versions.

* setuptools 7.0+
* pysam 0.7.5+
* numpy 1.6.2+
* scipy 0.9.0+
* docopt 0.6.1+
* cython 0.18+

Testing:
=========
Unit tests in the tests/ folder, run "nosetests test_rnacounter.py".

Testing files in the testfiles/ folder:
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after it.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.

Example::

    rnacounter testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf

The BAM contains 4041 reads all aligning perfectly on Gapdh (ENSMUSG00000057666) exons,
mostly on ENSMUSE00000487077 but also ENSMUSE00000751942 and ENSMUSE00000886744.
Nothing on other exons, which makes it a good example of badly conditioned input data...

The least squares method returns counts on the following transcripts:
ENSMUST00000117757, ENSMUST00000118875, ENSMUST00000147954
and nothing on ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588 .

Returns a count of 2459.62 (1091.71 RPK) for the gene.

