"""
Tool for generating star catalogues from the raw data collected by the
GALEX ultraviolet telescope.

For details of each available action use 'glcat <action> --help'.
"""

import argparse
import sys

from typing import NoReturn


def main() -> NoReturn:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.print_help()
    sys.exit(1)


if __name__ == "__main__":
    main()
