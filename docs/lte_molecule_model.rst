LTE Molecule Model
==================


There is a tool for modeling the full LTE spectrum of a molecule.
It uses the JPL or CDMS databases through their astroquery interfaces.

A very simple example looks like this:

.. code-block:: python
    
   import astropy.units as u
   from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_fitter, generate_model

   freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
                                                            fmin=90*u.GHz,
                                                            fmax=100*u.GHz)
   def modfunc(xarr, vcen, width, tex, column):
       return generate_model(xarr, vcen, width, tex, column, freqs=freqs, aij=aij,
                             deg=deg, EU=EU, partfunc=partfunc)

   fitter = generate_fitter(modfunc, name="CH3OH")

Details can be found in the API documentation:

.. automodule:: pyspeckit.spectrum.models.lte_molecule
    :members:
