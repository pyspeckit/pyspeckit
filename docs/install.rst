Installation and Requirements
=============================

.. hint::
    You can *pip_install* pyspeckit to get the latest release version: ::


        pip install pyspeckit



-------

PySpecKit requires at least the basic scientific packages:

* `numpy <https://numpy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net>`_
* `mpfit is included <https://github.com/segasai/astrolibpy/tree/master/mpfit>`_
* `scipy <https://www.scipy.org/>`_ is optional. It is  only required for RADEX
  grid interpolation and certain types of optimization
* python2.7 or `ordereddict <https://pypi.org/project/ordereddict>`_ for model parameter storage

You'll most likely want at least one of the following packages
to enable `file reading <readers>`_

* astropy_ >=0.4
* `atpy <http://atpy.github.com/>`_ (which depends on `asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_ [`github link <https://github.com/taldcroft/asciitable>`_] )
* `hdf5 <https://www.pytables.org/>`_

If you have pip (see https://pypi.org/project/pyspeckit), you can install with::

    pip install pyspeckit

Or the most recent version with either of these commands::

    pip install https://github.com/pyspeckit/pyspeckit/archive/master.zip

You can acquire the code with this clone command (see also :doc:`install_via_GitHub`)::

    git clone git@github.com:pyspeckit/pyspeckit.git pyspeckit
    cd pyspeckit
    python setup.py install

Or you can `Download the latest zip version <https://github.com/pyspeckit/pyspeckit/archive/refs/heads/master.zip>`_, 
then extract and install using the standard python method (but the pip install version of this is easier)::

    wget --no-check-certificate https://github.com/pyspeckit/pyspeckit/archive/refs/heads/master.zip
    unzip master.zip
    cd pyspeckit-pyspeckit-[commit]
    python setup.py install


You can also check out the `source code <https://github.com/pyspeckit/pyspeckit>`_

.. note ::
    If you use `easy_install pyspeckit` with the Enthought Python Distribution, you will
    most likely get a SandboxViolation error.  You can get around this by using `python
    setup.py install` or `pip install pyspeckit`.

.. toctree::
   install_via_GitHub
