Installation and Requirements
=============================

.. hint::
    You can *easy_install* or *pip_install* pyspeckit: ::


        pip install pyspeckit
        easy_install pyspeckit



-------

PySpecKit requires at least the basic scientific packages:

* `numpy <http://numpy.scipy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net>`_
* `mpfit is included <http://code.google.com/p/astrolibpy/source/browse/trunk/mpfit>`_
* `scipy <http://www.scipy.org/>`_ is optional. It is  only required for RADEX grid interpolation
* python2.7 or `ordereddict <http://pypi.python.org/pypi/ordereddict>`_ for model parameter storage

You'll most likely want at least one of the following packages
to enable `file reading <readers>`_

* `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits/Download>`_ >=2.4
* `atpy <http://atpy.github.com/>`_ (which depends on `asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_ [`github link <https://github.com/taldcroft/asciitable>`_] )
* `hdf5 <http://www.pytables.org/moin>`_

If you have pip (see http://pypi.python.org/pypi/pyspeckit), you can install with::

    pip install pyspeckit

Or the most recent version::

    pip install https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/get/tip.tar.gz

You can acquire the code with this clone command::

    hg clone https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org pyspeckit
    cd pyspeckit
    python setup.py install

or if you use git::

    git clone git://github.com/keflavich/pyspeckit.git pyspeckit
    cd pyspeckit
    python setup.py install

Or you can `Download the latest tarball version <https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/get/tip.tar.gz>`_, 
then extract and install using the standard python method::

    wget https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/get/tip.tar.gz
    tar -xzf pyspeckit-pyspeckit.bitbucket.org-tip.tar.gz
    cd pyspeckit-pyspeckit.bitbucket.org-tip
    python setup.py install


You can also check out the `source code <https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/src>`_

.. note ::
    If you use `easy_install pyspeckit` with the Enthought Python Distribution, you will
    most likely get a SandboxViolation error.  You can get around this by using `python
    setup.py install` or `pip install pyspeckit`.

