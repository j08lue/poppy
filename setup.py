#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
import glob

setup(name = 'poppy',
      version = '0.1.0',
      author = 'Jonas Bluethgen',
      author_email = 'bluthgen@nbi.ku.dk',
      packages = ['poppy'],
      url = 'http://www.gfy.ku.dk/~bluthgen',
      license = 'LICENSE.txt',
      description = 'Python tools for analyzing output data from the Parallel Ocean Program (POP).',
      long_description = open('README.md').read(),
      scripts = glob.glob(os.path.join('scripts', '*')),
      install_requires = [
          'oceanpy',
          'numpy',
          'matplotlib',
          'netCDF4',
          'pandas',
          ],
      dependency_links = [
          'https://github.com/j08lue/oceanpy/tarball/master#egg=oceanpy-0.1.0',
          ],
      )
