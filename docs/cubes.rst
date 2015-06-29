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

The most interesting features of the
:class:`~pyspeckit.cubes.SpectralCube.Cube` object are the
:meth:`~pyspeckit.cubes.SpectralCube.Cube.fiteach` method, which fits a model
spectrum to each element of the cube, and
:class:`mapplot <pyspeckit.cubes.mapplot.MapPlotter>`, which plots up various
projections of the cube.

:class:`Cube.mapplot <pyspeckit.cubes.mapplot.MapPlotter>` will create an interactive plot window.
You can click on any pixel shown in that window and pull up a second window
showing the spectrum at that pixel.  If you've fitted the cube, the associated
best-fit model will also be shown.
This interactive setup can be a bit fragile, though, so please report bugs
aggressively so we can weed them out!

The interactive viewer has a few button interactions described :meth:`here
<pyspeckit.cubes.mapplot.MapPlotter.mapplot>`.


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
