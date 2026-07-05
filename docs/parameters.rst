Parameters
==========

Model parameters are very flexible, and can be accessed and modified in many
parallel ways.

The `parinfo` class is built on top of `lmfit-py's parameters
<https://github.com/lmfit/lmfit-py/blob/master/doc/parameters.rst>`_ for compatibility
with lmfit-py, but it builds on that.  The code for the parameter overloading is in
`parinfo.py <https://github.com/pyspeckit/pyspeckit/blob/master/pyspeckit/spectrum/parinfo.py>`_.

.. _constraining-parameters:

Constraining fit parameters: limits, fixed, and tied
----------------------------------------------------

Every call to ``sp.specfit(...)`` accepts, in addition to ``guesses``, a set
of per-parameter constraint keywords.  Each is a list with one entry per
parameter (i.e., ``npars * npeaks`` entries; for a two-Gaussian fit, six
entries in the order amplitude, center, width, amplitude, center, width):

``limits``
    A list of two-element ``(min, max)`` tuples specifying the lower and
    upper bounds on each parameter.

``limited``
    A list of two-element ``(bool, bool)`` tuples that turn the
    corresponding entries of ``limits`` on or off.  **Entries in ``limits``
    are ignored unless the matching ``limited`` flag is `True`**; the common
    combination ``limits=(0,0), limited=(False,False)`` therefore means
    "unbounded", not "pinned to zero".

``fixed``
    A list of booleans.  A `True` entry holds that parameter exactly at its
    input guess during the fit.

``tied``
    A list of strings.  An empty string means "free"; a non-empty string is
    an expression evaluated in terms of the other parameters, referred to as
    ``p[0]``, ``p[1]``, ... (e.g. ``'p[1]+14.38'`` forces this parameter to
    equal parameter 1 plus a constant offset).

These work with both the default `mpfit` back-end and with lmfit-py
(``use_lmfit=True``).

Setting parameter limits (boundary conditions)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following is a complete, runnable example.  The synthetic line has a true
amplitude of 5, but the fit constrains the amplitude to the range ``[0, 3]``,
so the fitted amplitude pegs at exactly 3.0:

.. code-block:: python

    import numpy as np
    import pyspeckit

    xaxis = np.linspace(-50., 150., 100)
    sigma = 10.
    center = 50.
    amplitude = 5.

    np.random.seed(0)
    synth_data = amplitude * np.exp(-(xaxis - center)**2 / (2 * sigma**2))
    stddev = 0.1
    noise = np.random.randn(xaxis.size) * stddev
    error = stddev * np.ones_like(synth_data)
    data = noise + synth_data

    sp = pyspeckit.Spectrum(data=data, error=error, xarr=xaxis,
                            xarrkwargs={'unit': 'km/s'},
                            unit='K')

    # define shorthand first:
    T, F = True, False
    sp.specfit(fittype='gaussian', guesses=[1, 45, 5],
               limited=[(T, T), (F, F), (T, F)],
               limits=[(0, 3), (0, 0), (0, 0)])
    print(sp.specfit.parinfo)

The amplitude is bounded on both sides (``limited=(T,T)``) to the range 0-3,
the center is unconstrained, and the width has only a lower bound of zero
(``limited=(T,F)``, so the second element of its ``limits`` entry is
ignored).  The output shows the amplitude clamped at the upper boundary::

    Param #0   AMPLITUDE0 =            3 +/-               0   Range:     [0,3]
    Param #1       SHIFT0 =      50.0044 +/-        0.183198
    Param #2       WIDTH0 =      13.2506 +/-        0.149581   Range:   [0,inf)

(A parameter that ends up pegged at a boundary has zero derivative there, so
its reported error is 0 - a useful hint that the limit is active.)

One-sided limits are useful for forcing a component to be an emission or an
absorption line:

.. code-block:: python

    sp.specfit(fittype='gaussian', guesses=[-1,5,1,0.5,2,1],
        limits=[(0,0), (0,0), (0,0), (0,0), (0,0), (0,0)],
        limited=[(F,T), (F,F), (T,F), (T,F),  (F,F), (T,F)])

In this example, there are two gaussian components being fitted because a
Gaussian takes 3 parameters, an amplitude, a center, and a width, and there
are 6 parameters in the input guesses.

The first line is forced to be an absorption line: its limits are ``(0,0)``
but ``limited=(F,T)`` so only the 2nd element, the upper limit, is respected:
the amplitude is forced to be :math:`A\leq 0`.

The second line has its amplitude (the 4th parameter in ``guesses``) forced
*positive* since its limits are also ``(0,0)`` but its ``limited=(T,F)``.

Both lines have their *widths* forced to be positive, which is true by
default: there is no meaning to a negative width, since the width enters into
the equation for a gaussian as :math:`\sigma^2`.

Note that the need to limit parameters is the main reason for the existence
of `lmfit-py <https://github.com/lmfit/lmfit-py>`_ and `mpfit
<https://github.com/segasai/astrolibpy/blob/4993aa4e7c1001efe7c00048ec2b9d5ccac83ff7/mpfit/mpfit.py>`_.

Fixing parameters
^^^^^^^^^^^^^^^^^

To hold a parameter exactly at its input guess, use ``fixed``.  Continuing
the example above, this fits only the amplitude and center while the width is
held at 8:

.. code-block:: python

    sp.specfit(fittype='gaussian', guesses=[4, 45, 8],
               fixed=[F, F, T])
    print(sp.specfit.parinfo)

::

    Param #0   AMPLITUDE0 =      5.46492 +/-       0.0377455
    Param #1       SHIFT0 =      50.0387 +/-       0.0781422
    Param #2       WIDTH0 =            8 (fixed)  Range:   [0,inf)

Tying parameters
^^^^^^^^^^^^^^^^

``tied`` expresses one parameter as a function of the others.  For example,
to fit the [S II] 6716,6731 doublet with the line separation locked to its
laboratory value:

.. code-block:: python

    xaxis2 = np.linspace(6690., 6760., 200)
    sigma2 = 1.5
    synth_data2 = (1.5 * np.exp(-(xaxis2 - 6716.44)**2 / (2 * sigma2**2)) +
                   1.0 * np.exp(-(xaxis2 - 6730.82)**2 / (2 * sigma2**2)))
    np.random.seed(1)
    noise2 = np.random.randn(xaxis2.size) * 0.05
    error2 = 0.05 * np.ones_like(synth_data2)
    sp2 = pyspeckit.Spectrum(data=synth_data2 + noise2, error=error2,
                             xarr=xaxis2,
                             xarrkwargs={'unit': 'angstrom'},
                             unit='erg/s/cm2/AA')

    sp2.specfit(fittype='gaussian', guesses=[1, 6716, 1, 1, 6731, 1],
                tied=['', '', '', '', 'p[1]+14.38', ''])
    print(sp2.specfit.parinfo)

::

    Param #0   AMPLITUDE0 =      1.48456 +/-       0.0220743
    Param #1       SHIFT0 =      6716.46 +/-       0.0218963
    Param #2       WIDTH0 =      1.52732 +/-       0.0262233   Range:   [0,inf)
    Param #3   AMPLITUDE1 =     0.980497 +/-       0.0220252
    Param #4       SHIFT1 =      6730.84 +/-               0  Tied: p[1]+14.38
    Param #5       WIDTH1 =      1.53413 +/-       0.0397928   Range:   [0,inf)

The second component's center is exactly 14.38 Angstroms redward of the
first, by construction.  Tied parameters can be combined freely with
``limits``/``limited`` on the other parameters.  :doc:`example_sdss` shows a
more complete example using ``tied``.

If you use `lmfit-py <https://github.com/lmfit/lmfit-py>`_ by specifying
``use_lmfit=True``, you can use the more advanced `mathematical constraints
<https://lmfit.github.io/lmfit-py/constraints.html>`_ permitted by lmfit-py.

Inspecting and reusing the ``parinfo``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After any fit, ``sp.specfit.parinfo`` is a
:class:`~pyspeckit.spectrum.parinfo.ParinfoList` holding the fitted values,
uncertainties, and all of the constraint metadata.  You can index it by
parameter name or number:

.. code-block:: python

    sp.specfit.parinfo               # pretty-printed table of all parameters
    sp.specfit.parinfo['WIDTH0']     # a single Parinfo object
    sp.specfit.parinfo.values        # [5.46, 50.04, 8.0]
    sp.specfit.parinfo.errors        # [0.038, 0.078, 0.0]

    pi = sp.specfit.parinfo
    pi['AMPLITUDE0'].value, pi['AMPLITUDE0'].error, pi['AMPLITUDE0'].limits

You can also modify a ``parinfo`` and pass it back in as the starting point
(and constraint specification) for a new fit, instead of using the keyword
lists:

.. code-block:: python

    pi = sp.specfit.parinfo
    pi['AMPLITUDE0'].limited = (True, True)
    pi['AMPLITUDE0'].limits = (0, 3)
    pi['AMPLITUDE0'].value = 2.5
    pi['WIDTH0'].fixed = False
    sp.specfit(parinfo=pi)

Model-specific keywords: ``minpars``, ``maxpars``, ``limitedmin``, ``limitedmax``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As an alternative to the tuple-based ``limits``/``limited`` keywords,
``specfit`` also accepts the split form used by several of the model
wrappers: ``minpars`` and ``maxpars`` are lists of scalar lower/upper
bounds, and ``limitedmin``/``limitedmax`` are the corresponding lists of
booleans.  ``minpars=[0,0,0], limitedmin=[True,False,True]`` is equivalent to
``limits=[(0,0),(0,0),(0,0)], limited=[(True,False),(False,False),(True,False)]``.

Some models apply sensible physical limits by default.  The ammonia model,
for example (parameters ``tkin, tex, ntot, width, xoff_v, fortho``), uses
default constraints equivalent to passing

.. code-block:: python

    limitedmin=(True, True, True, True, False, True),
    limitedmax=(False, False, True, False, False, True),
    minpars=(2.7315, 2.7315, 5, 0, 0, 0),
    maxpars=(0, 0, 25, 0, 0, 1)

i.e., temperatures bounded below by the CMB temperature, log column density
between 5 and 25, positive linewidth, and ortho fraction between 0 and 1.
You can override any of these by passing your own values for those keywords
to ``specfit``.

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

