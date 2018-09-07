LTE Molecule Model
==================


There is a tool for modeling the full LTE spectrum of a molecule.
It uses the CDMS / VAMDC database
(http://portal.vamdc.eu/vamdc_portal/home.seam) and the vamdclib
(http://vamdclib.readthedocs.io) python library to compute the
partition function of a molecule.  It uses astroquery.splatalogue
(http://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html)
to identify the line frequencies and energy levels.

A very simple example looks like this:

.. code-block:: python


   freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
                                                            fmin=90*u.GHz,
                                                            fmax=100*u.GHz)
   def modfunc(xarr, vcen, width, tex, column):
       return generate_model(xarr, vcen, width, tex, column, freqs=freqs, aij=aij,
                             deg=deg, EU=EU, partfunc=partfunc)

   fitter = generate_fitter(modfunc, name="CH3OH")

The molecular parameter lookup stage is often slow and may be a bottleneck.

Details can be found in the API documentation:

.. automodule:: pyspeckit.spectrum.models.lte_molecule
    :members:
