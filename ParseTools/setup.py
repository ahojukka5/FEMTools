# Copyright (c) 2012 FILL IN
from distutils.core import setup
#from setuptools import setup

from parsetoolslib import meta

setup(name='ParseTools',
      version=meta.__version__,
      author=meta.__author__,
      description='FILL IN',
      scripts=['bin/parsetools'],
      package_dir={'parsetoolslib':'parsetoolslib'},
      packages=['parsetoolslib'],
)
