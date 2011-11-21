#!/usr/bin/env python

import os, shutil
from distutils.core import setup

setup(name='pyspeckit',
      version='0.1',
      description='Spectroscopic Toolkit for fitting and manipulating spectroscopic data in python',
      author=['Adam Ginsburg','Jordan Mirocha'],
      author_email=['adam.g.ginsburg@gmail.com','mirochaj@gmail.com'],
      url='https://bitbucket.org/keflavich/spectroscopic-toolkit-astronomy/',
      packages=['pyspeckit','pyspeckit.spectrum','pyspeckit.spectrum.models','pyspeckit.spectrum.readers','pyspeckit.spectrum.writers','pyspeckit.spectrum.speclines',
          'pyspeckit.cubes', 'pyspeckit.wrappers','mpfit'],
      package_dir={'pyspeckit.spectrum.speclines':'pyspeckit/spectrum/speclines','mpfit':'mpfit'},
      package_data={'pyspeckit.spectrum.speclines':['splatalogue.csv']},
     )

# Copy default config file to directory $HOME/.pyspeckit
HOME = os.environ.get('HOME')
CWD = os.getcwd()
if not os.path.exists('%s/.pyspeckit' % HOME):
    os.mkdir('%s/.pyspeckit' % HOME)
    
if not os.path.exists('%s/.pyspeckit/config' % HOME):
    shutil.copyfile('%s/pyspeckit/config_default' % CWD, '%s/.pyspeckit/config' % HOME)
