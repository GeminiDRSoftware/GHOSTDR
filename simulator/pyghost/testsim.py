import pyghost.tests
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-cosmics', action='store_false', help="*don't* add in"
                        " cosmic rays to the generated frames")
    parser.add_argument('-crplane', action='store_true', help='output cosmic'
                        ' ray locations as separate fits files')
    parser.add_argument('-hpplane', action='store_true', help='output hot pixel'
                        ' locations as separate fits files')
    parser.add_argument('-split', action='store_true', help='output frames as'
                        ' individual files instead of bundled together (as a'
                        ' MEF)')
    parser.add_argument('-check', action='store_true', help='output additional'
                        ' check files that contain the input spectrum, interpolated'
                        ' onto the pixel grid for each order')
    # will raise if too few/many args or unsupported options used
    args = parser.parse_args()
    pyghost.tests.run(**vars(args))  # pylint: disable=star-args
