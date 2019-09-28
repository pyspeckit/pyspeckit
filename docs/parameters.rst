Parameters
==========

Model parameters are very flexible, and can be accessed and modified in many
parallel ways.

The `parinfo` class is built on top of `lmfit-py's parameters
<https://github.com/lmfit/lmfit-py/blob/master/doc/parameters.rst>`_ for compatibility
with lmfit-py, but it builds on that.  The code for the parameter overloading is in
`parinfo.py <https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/src/tip/pyspeckit/spectrum/parinfo.py>`_.

Simple Example
--------------
Start with a simple example: if you want to limit parameters to be within some range, use
the `limits` and `limited` parameters.

.. code-block:: python

    # define shorthand first:
    T,F = True,False
    sp.specfit(fittype='gaussian', guesses=[-1,5,1,0.5,2,1],
        limits=[(0,0), (0,0), (0,0), (0,0), (0,0), (0,0)],
        limited=[(F,T), (F,F), (T,F), (T,F),  (F,F), (T,F)])

In this example, there are two gaussian components being fitted because a
Gaussian takes 3 parameters, an amplitude, a center, and a width, and there are
6 parameters in the input guesses.

The first line is forced to be an absorption line: its limits are `(0,0)` but
`limited=(F,T)` so only the 2nd parameter, the upper limit, is respected: the amplitude
is forced to be :math:`A\leq 0`.  

The second line has its amplitude (the 4th parameter in `guesses`) forced
*positive* since its limits are also `(0,0)` but its `limited=(T,F)`.  

Both lines have their *widths* forced to be positive, which is true by default:
there is no meaning to a negative width, since the width enters into the
equation for a gaussian as :math:`\sigma^2`.

Note that the need to limit parameters is the main reason for the existence of `lmfit-py <https://github.com/lmfit/lmfit-py>`_
and `mpfit <https://github.com/segasai/astrolibpy/blob/4993aa4e7c1001efe7c00048ec2b9d5ccac83ff7/mpfit/mpfit.py>`_. 

Tying Parameters
----------------
It is also possible to explicitly state that one parameter depends on another.
If, for example, you want to fit two gaussians, but they must be at a fixed
wavelength separation from one another (e.g., for fitting the [S II] doublet),
use `tied`:

.. code-block:: python

    sp.specfit(fittype='gaussian', guesses=[1,6716,1,0.5,6731,1],
        tied=['','','','','p[1]',''])

If you use `lmfit-py` by specifying `use_lmfit=True`, you can use the more advanced `mathematical constraints
<http://lmfit.github.com/lmfit-py/constraints.html>`_ permitted by lmfit-py.

:doc:`example_sdss` shows a more complete example using `tied`.

Making your own `parinfo`
-------------------------
You can also build a `parinfo` class directly.
Currently, the best example of this is in `tests/test_formaldehyde_mm_radex.py`.

Here's an example of how you would set up a fit using `parinfo` directly.

.. WARNING:: 

    There is a bug in the use_lmfit section of this code that keeps it from working properly.  =(

.. code-block:: python

    amplitude0 = pyspeckit.parinfo.Parinfo(n=0, parname='Amplitude0',
        shortparname='$A_0$', value=1, limits=[0, 100], limited=(True,True)) 
    width0 = pyspeckit.parinfo.Parinfo(n=2, parname='Width0',
        shortparname='$\sigma_0$', value=1, limits=(0, 0), limited=(True,False))
    center0 = pyspeckit.parinfo.Parinfo(n=1, parname='Center0',
        shortparname='$\Delta x_0$', value=6716, limits=(0, 0), limited=(False,False))
    amplitude1 = pyspeckit.parinfo.Parinfo(n=3, parname='Amplitude1',
        shortparname='$A_1$', value=1, limits=[0, 100], limited=(True,True)) 
    width1 = pyspeckit.parinfo.Parinfo(n=5, parname='Width1',
        shortparname='$\sigma_1$', value=1, limits=(0, 0), limited=(True,False))
    center1 = pyspeckit.parinfo.Parinfo(n=4, parname='Center1',
        shortparname='$\Delta x_1$', value=6731, limits=(0, 0),
        limited=(False,False), tied=center0)

    parinfo = pyspeckit.parinfo.ParinfoList([amplitude0,center0,width0,amplitude1,center1,width1])

    sp.specfit(parinfo=parinfo, use_lmfit=True)

