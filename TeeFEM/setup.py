# Copyright (c) 2012 Jukka Aho

from distutils.core import setup
from teefem import meta

#from setuptools import setup

setup(name='TeeFEM',
      version=meta.__version__,
      author=meta.__author__,
	  author_email = meta.__email__,
      description='TeeFEM - Teekkarin FEM',
      scripts=['bin/teefem'],
      #package_dir={'teefemlib':'teefemlib'},
      packages=['teefem','teefem.models'],
)
