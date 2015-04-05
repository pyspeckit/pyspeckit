Cubes
-----

Pyspeckit can do a few things with spectral cubes.  The most interesting is the
spectral line fitting. 

:class:`~pyspeckit.cubes.SpectralCube.Cube` objects have a
:meth:`~pyspeckit.cubes.SpectralCube.Cube.fiteach` method that will fit each
spectral line within a cube.  It can be made to do this in parallel with the
``multicore`` option.

As of version 0.16, pyspeckit cubes can be read from SpectralCube_ objects:

.. code::

   >>> pcube = pyspeckit.Cube(cube=mySpectralCube)

Otherwise, they can be created from FITS cubes on disk:

.. code::

   >>> pcube = pyspeckit.Cube(filename="mycube.fits")

or from arrays:

.. code::

   >>> mycube = np.random.randn(250,50,50)
   >>> myxaxis = np.linspace(-100,100,250)
   >>> pcube = pyspeckit.Cube(cube=mycube, xarr=myxaxis, xunit='km/s')


.. automodule:: pyspeckit.cubes.SpectralCube
.. autoclass:: Cube
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:
.. autoclass:: CubeStack
    :show-inheritance:
    :members:


.. automodule:: pyspeckit.cubes.mapplot
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:

.. automodule:: pyspeckit.cubes.cubes
    :members:

.. _SpectralCube: http://spectral-cube.readthedocs.org/en/latest/
