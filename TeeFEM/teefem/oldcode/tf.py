# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 19:43:07 2012

Jukka Aho <jukka.aho@kapsi.fi>

TeeFEM -- TeekkariFEM

Ratkaisee tasojännitystilan probleemia.

http://assets.en.oreilly.com/1/event/27/Best%20practices%20for%20_scripting_%20with%20Python%203%20Paper.pdf

"""

import sys
import optparse
import meta
import logging
logging.basicConfig(level = logging.DEBUG, filename='teefem.log')
logging.debug('Starting log')

def test(*args, **kwds):
    print("Running tests")
    import unittest
    import tests
    unittest.main()

def main(prog_args):
    
    functions = {'test': test}    
    
    parser = optparse.OptionParser(version=meta.__version__)
    opt, args = parser.parse_args(prog_args)
    print opt,args
    for arg in args:
        if arg in functions:
            functions[arg]()

if __name__ == '__main__':
    test()
#    sys.exit(main(sys.argv))

