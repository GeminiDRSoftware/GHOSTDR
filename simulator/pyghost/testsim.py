import pyghost.tests
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-crplane', action='store_true', help='output cosmic '
                        'ray locations as separate fits files')
    parser.add_argument('-split', action='store_true', help='output frames as '
                        'individual files instead of bundled together (as a '
                        'MEF)')
    # will raise if too few/many args or unsupported options used
    args = parser.parse_args()
    pyghost.tests.run(crplane=args.crplane, split=args.split)
