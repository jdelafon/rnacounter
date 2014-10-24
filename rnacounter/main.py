
try:
    from rnacounter import rnacounter
except ImportError:
    import rnacounter
import docopt

def main():
    args = docopt.docopt(rnacounter.usage_string())
    if args['join']:
        rnacounter.join([args['TAB']]+args['TAB2'])
    else:
        bamname, annotname, options = rnacounter.parse_args(args)
        rnacounter.rnacounter_main(bamname,annotname, options)

if __name__ == '__main__':
    main()
