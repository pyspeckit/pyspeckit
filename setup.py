#!/usr/bin/env python

from distutils.core import setup

setup(name='pyspeckit',
      version='0.1',
      description='Spectroscopic Toolkit for fitting and manipulating spectroscopic data in python',
      author=['Adam Ginsburg','Jordan Mirocha'],
      author_email=['adam.g.ginsburg@gmail.com','mirochaj@gmail.com'],
      url='https://bitbucket.org/keflavich/spectroscopic-toolkit-astronomy/',
      packages=['pyspeckit','pyspeckit.spectrum','pyspeckit.spectrum.models','pyspeckit.spectrum.readers','pyspeckit.spectrum.writers','pyspeckit.spectrum.speclines',
          'pyspeckit.cubes', 'pyspeckit.wrappers'],
      package_dir={'spectrum.speclines':'spectrum/speclines'},
      package_data={'spectrum.speclines':['splatalogue.csv']},
     )
