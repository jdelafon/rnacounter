
try:
    from rnacounter import rnacounter
except ImportError:
    import rnacounter
import docopt


# Copy the version given in setup.py
import pkg_resources
__version__ = pkg_resources.require("rnacounter")[0].version


def main():
    args = docopt.docopt(rnacounter.usage_string(), version=__version__)
    if args['join']:
        rnacounter.join([args['TAB']]+args['TAB2'])
    else:
        bamname, annotname, options = rnacounter.parse_args(args)
        rnacounter.rnacounter_main(bamname,annotname, options)

if __name__ == '__main__':
    main()
