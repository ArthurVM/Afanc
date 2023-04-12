#! /usr/bin/env python3
"""
Command line interface for Afanc metagenomics disambiguator.
"""
import sys

from Afanc.parser import base_parser, initLogFiles


def main():
    args = base_parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        base_parser.print_help()


if __name__ == "__main__":
    main()
