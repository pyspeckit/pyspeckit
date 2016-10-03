Installation and Requirements
=============================

.. hint::
    You can *easy_install* or *pip_install* pyspeckit to get the latest release version: ::


        pip install pyspeckit
        easy_install pyspeckit



-------

PySpecKit requires at least the basic scientific packages:

* `numpy <http://numpy.scipy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net>`_
* `mpfit is included <https://github.com/segasai/astrolibpy/tree/master/mpfit>`_
* `scipy <http://www.scipy.org/>`_ is optional. It is  only required for RADEX
  grid interpolation and certain types of optimization
* python2.7 or `ordereddict <http://pypi.python.org/pypi/ordereddict>`_ for model parameter storage

You'll most likely want at least one of the following packages
to enable `file reading <readers>`_

* astropy_ >=0.4
* `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits/Download>`_ >=2.4
* `atpy <http://atpy.github.com/>`_ (which depends on `asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_ [`github link <https://github.com/taldcroft/asciitable>`_] )
* `hdf5 <http://www.pytables.org/moin>`_

If you have pip (see http://pypi.python.org/pypi/pyspeckit), you can install with::

    pip install pyspeckit

Or the most recent version::

    pip install https://bitbucket.org/pyspeckit/pyspeckit/get/master.tar.gz

You can acquire the code with this clone command (see also :doc:`install_via_GitHub`)::

    git clone git@bitbucket.org:pyspeckit/pyspeckit.git pyspeckit
    cd pyspeckit
    python setup.py install

Or you can `Download the latest tarball version <https://bitbucket.org/pyspeckit/pyspeckit/get/master.tar.gz>`_, 
then extract and install using the standard python method (but the pip install version of this is easier)::

    wget --no-check-certificate https://bitbucket.org/pyspeckit/pyspeckit/get/master.tar.gz
    tar -xzf master.tar.gz
    cd pyspeckit-pyspeckit-[commit]
    python setup.py install


You can also check out the `source code <https://bitbucket.org/pyspeckit/pyspeckit/src>`_

.. note ::
    If you use `easy_install pyspeckit` with the Enthought Python Distribution, you will
    most likely get a SandboxViolation error.  You can get around this by using `python
    setup.py install` or `pip install pyspeckit`.

.. note ::
   pyspeckit is hosted on both bitbucket and github.  Both versions are kept up
   to date, so it should not matter which one you choose to install.
