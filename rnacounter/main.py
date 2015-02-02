#!/usr/bin/env python

# Entry point

try:
    from rnacounter import rnacounter
except ImportError:
    import rnacounter
import docopt, os
from pkg_resources import resource_filename, require


# Copy the version given in setup.py
__version__ = require("rnacounter")[0].version

def main():
    args = docopt.docopt(rnacounter.usage_string(), version=__version__)
    if args['join']:
        rnacounter.join([args['TAB']]+args['TAB2'], args.get('--output'))
    elif args['test']:
        options = rnacounter.parse_args(args)
        bamname = os.path.abspath(resource_filename('testfiles', 'gapdhKO.bam'))
        annotname = os.path.abspath(resource_filename('testfiles', 'mm9_3genes_renamed.gtf'))
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
