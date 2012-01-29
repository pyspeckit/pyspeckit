#!/usr/bin/env python

import os, shutil
from distutils.core import setup

with open('README.txt') as file:
    long_description = file.read()

version_base="0.1"

if os.path.exists(".hg"):
    try:
        import subprocess
        currentversion = subprocess.Popen(["hg","id","--num"],stdout=subprocess.PIPE).communicate()[0].strip().strip("+")
        version = version_base+"hg"+currentversion
    except:
        # is this bad practice?  I don't care if it's an import error, attribute error, or value error...
        version = version_base

setup(name='pyspeckit',
      version=version,
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
      requires=['matplotib','numpy'],
      classifiers=[
                   "Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                  ],
      
     )

# Copy default config file to directory $HOME/.pyspeckit
HOME = os.environ.get('HOME')
CWD = os.getcwd()
if not os.path.exists('%s/.pyspeckit' % HOME):
    os.mkdir('%s/.pyspeckit' % HOME)
    
if not os.path.exists('%s/.pyspeckit/config' % HOME):
    shutil.copyfile('%s/pyspeckit/config_default' % CWD, '%s/.pyspeckit/config' % HOME)
