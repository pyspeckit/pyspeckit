#!/usr/bin/env python

import os, shutil
from distutils.core import setup

with open('README.txt') as file:
    long_description = file.read()

setup(name='pyspeckit',
      version='0.1',
      description='Toolkit for fitting and manipulating spectroscopic data in python',
      long_description=long_description,
      author=['Adam Ginsburg','Jordan Mirocha'],
      author_email=['adam.g.ginsburg@gmail.com', 'mirochaj@gmail.com',
          'pyspeckit@gmail.com'], 
      url='https://pyspeckit.bitbucket.org/',
      packages=['pyspeckit', 'pyspeckit.spectrum', 'pyspeckit.spectrum.models',
          'pyspeckit.spectrum.readers', 'pyspeckit.spectrum.writers',
          'pyspeckit.spectrum.speclines', 'pyspeckit.cubes',
          'pyspeckit.wrappers','mpfit','parallel_map'],
      package_dir={'pyspeckit.spectrum.speclines':'pyspeckit/spectrum/speclines',
          'mpfit':'mpfit'}, 
      package_data={'pyspeckit.spectrum.speclines':['splatalogue.csv']},
     )

# Copy default config file to directory $HOME/.pyspeckit
HOME = os.environ.get('HOME')
CWD = os.getcwd()
if not os.path.exists('%s/.pyspeckit' % HOME):
    os.mkdir('%s/.pyspeckit' % HOME)
    
if not os.path.exists('%s/.pyspeckit/config' % HOME):
    shutil.copyfile('%s/pyspeckit/config_default' % CWD, '%s/.pyspeckit/config' % HOME)
