#!/usr/bin/env python
# Copyright (c) 2012 FILL IN

import sys
import optparse
import meta

from jukka import *

def main(prog_args):
    parser = optparse.OptionParser(version=meta.__version__)
    opt, args = parser.parse_args(prog_args)

if __name__ == '__main__':
    sys.exit(main(sys.argv))

