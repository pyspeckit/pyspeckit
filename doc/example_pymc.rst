
.. include:: <isogrk3.txt>

Monte Carlo examples
~~~~~~~~~~~~~~~~~~~~
There are (at least) two packages implementing Monte Carlo sampling available
in python: `pymc <http://code.google.com/p/pymc/>`_ and `emcee
<http://danfm.ca/emcee/>`_.  `pyspeckit` implements interfaces to both.  With
the pymc interface, it is possible to define priors that strictly limit the
parameter space.  So far that is not possible with emcee.

The examples below use a custom plotting package from `agpy
<http://code.google.com/p/agpy/>`_.  It is a relatively simple but convenient
wrapper around numpy's histogram2d.  `pymc_plotting
<http://code.google.com/p/agpy/source/browse/trunk/agpy/pymc_plotting.py>`_
takes care of indexing, percentile determination, and coloring.

:: 

    import pyspeckit

    # Create our own gaussian cetered at 0 with width 1, amplitude 5, and
    # gaussian noise with amplitude 1
    x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
    e = np.random.randn(50)
    d = np.exp(-np.asarray(x)**2/2.)*5 + e

    # create the spectrum object
    sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
    # fit it
    sp.specfit(fittype='gaussian')
    # then get the pymc values
    MCuninformed = sp.specfit.get_pymc()
    MCwithpriors = sp.specfit.get_pymc(use_fitted_values=True)
    MCuninformed.sample(5000,burn=1000,tune_interval=250)
    MCwithpriors.sample(5000,burn=1000,tune_interval=250)

    # MC vs least squares:
    print sp.specfit.parinfo
    # Param #0       HEIGHT =       0.0109281 +/-        0.195809 
    # Param #1   AMPLITUDE0 =         4.57169 +/-        0.651568 
    # Param #2       SHIFT0 =     0.000182704 +/-        0.205333 
    # Param #3       WIDTH0 =         1.23275 +/-        0.219416   Range:   [0,inf)
    

    print MCuninformed.stats()['AMPLITUDE0'],MCuninformed.stats()['WIDTH0']
    # {'95% HPD interval': array([ 2.80477894,  5.19230514]),
    #  'mc error': 0.0096630856947778534,
    #  'mean': 3.8550042004763823,
    #  'n': 4000,
    #  'quantiles': {2.5: 2.7001147437515614,
    #   25: 3.4136653207924721,
    #   50: 3.8216640136271343,
    #   75: 4.2658237064708775,
    #   97.5: 5.1241797879905224},
    #  'standard deviation': 0.61114720041777293}
    # {'95% HPD interval': array([ 0.72443531,  1.56684083]),
    #  'mc error': 0.0035649461062602994,
    #  'mean': 1.1868052003204772,
    #  'n': 4000,
    #  'quantiles': {2.5: 0.7886369815120019,
    #   25: 1.030731243691529,
    #   50: 1.1608275996051221,
    #   75: 1.3238167303614339,
    #   97.5: 1.6719731801069622},
    #  'standard deviation': 0.22546698863062389}
    
    print MCwithpriors.stats()['AMPLITUDE0'],MCwithpriors.stats()['WIDTH0']
    # {'95% HPD interval': array([ 3.14723955,  4.75558755]),
    #  'mc error': 0.0066693160876581236,
    #  'mean': 4.0120510164373284,
    #  'n': 4000,
    #  'quantiles': {2.5: 3.1905811738800782,
    #   25: 3.7258959755336791,
    #   50: 4.0192907032241969,
    #   75: 4.2848569265514413,
    #   97.5: 4.8263499154520932},
    #  'standard deviation': 0.42180458545205723},
    # {'95% HPD interval': array([ 0.84659241,  1.37188539]),
    #  'mc error': 0.0021143746687207888,
    #  'mean': 1.132295044520145,
    #  'n': 4000,
    #  'quantiles': {2.5: 0.87023370176467096,
    #   25: 1.0471063188439775,
    #   50: 1.1319776513946447,
    #   75: 1.2149336329380425,
    #   97.5: 1.4028230972640638},
    #  'standard deviation': 0.13372479560243336}


    # optional
    import agpy.pymc_plotting
    import pylab
    agpy.pymc_plotting.hist2d(MCuninformed,'AMPLITUDE0','WIDTH0',clear=True)
    agpy.pymc_plotting.hist2d(MCwithpriors,'AMPLITUDE0','WIDTH0',contourcmd=pylab.contour,colors=[(0,1,0,1),(0,0.5,0.5,1),(0,0.75,0.75,1),(0,1,1,1),(0,0,1,1)],clear=False)


    # Now do the same with emcee
    emcee_ensemble = sp.specfit.get_emcee()
    p0 = emcee_ensemble.p0 * (np.random.randn(*emcee_ensemble.p0.shape) / 10. + 1.0)
    pos,logprob,state = emcee_ensemble.run_mcmc(p0,4000)

    plotdict = {'AMPLITUDE0':emcee_ensemble.chain[:,:,1].ravel(),
                'WIDTH0':emcee_ensemble.chain[:,:,3].ravel()}
    agpy.pymc_plotting.hist2d(plotdict,'AMPLITUDE0','WIDTH0',fignum=2)

.. figure:: ../images/pymc_params_example.png 
    :alt: Examples of the pymc Monte Carlo two-dimensional parameter histogram
        (marginalized over X-offset and Y-offset) with and without priors; with
        priors is shown in contour lines
    :figwidth: 800
    :width: 800

.. figure:: ../images/emcee_params_example.png 
    :alt: Same as above, but using emcee
    :figwidth: 800
    :width: 800


