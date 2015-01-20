#!/usr/bin/env python

# Entry point

try:
    from rnacounter import rnacounter
except ImportError:
    import rnacounter
import docopt, os


# Copy the version given in setup.py
import pkg_resources
__version__ = pkg_resources.require("rnacounter")[0].version


def main():
    args = docopt.docopt(rnacounter.usage_string(), version=__version__)
    if args['join']:
        rnacounter.join([args['TAB']]+args['TAB2'])
    elif args['test']:
        options = rnacounter.parse_args(args)
        local = os.path.abspath(os.path.dirname(__file__))
        testfiles = os.path.join(local, "..", "tests", "testfiles")
        bamname = os.path.join(testfiles, "gapdhKO.bam")
        annotname = os.path.join(testfiles, "mm9_3genes_renamed.gtf")
        rnacounter.rnacounter_main(bamname,annotname, options)
    else:
        bamname = os.path.abspath(args['BAM'])
        annotname = os.path.abspath(args['GTF'])
        assert os.path.exists(bamname), "BAM file not found: %s" %bamname
        assert os.path.exists(annotname), "GTF file not found: %s" %annotname
        options = rnacounter.parse_args(args)
        rnacounter.rnacounter_main(bamname,annotname, options)


if __name__ == '__main__':
    main()
