Examples
--------

Check out the `flickr gallery <http://flic.kr/s/aHsjykNJWt>`_.

.. .. flickr:: 
..    :flickrID: 109615

Want your image or example included?  `E-mail us <mailto:pyspeckit@gmail.com>`_.

..
    Radio
    -----
     * :doc:`example_h2co`
     * :doc:`example_nh3`
     * :doc:`example_hcop`
    Near IR
    -------
     * :doc:`example_NIR_cube`
    Optical
    -------
     * :doc:`example_sdss`
    UV
    --
     * :doc:`fit_custom_model`
     


.. toctree::
   :maxdepth: 3

   Build a spectrum from scratch (and fit a gaussian) <example_fromscratch>
   Build a spectrum from scratch (and fit the continuum) <example_continuum_fromscratch>
   Minimal gaussian cube fitting <example_minimal_cube>
   Build a cube from scratch and fit many gaussians <example_cube_fromscratch>
   Radio H2CO <example_h2co>
   Radio H2CO mm lines <example_h2co_mm>
   Radio NH3: shows hyperfine fitting with fixed amplitude, width <example_nh3>
   Radio NH3 cube <example_nh3_cube>
   Radio NH3 - fit N hyperfines with no temperature <nh3_tau>
   Radio N2H+ <example_n2hp>
   Radio N2H+ cube <example_n2hp_cube>
   Radio ortho-NH2D <example_o-nh2d>
   Radio para-NH2D <example_p-nh2d>
   Radio HCN: shows hyperfine with varying amplitude & width <example_hcn>
   Radio HCO+ <example_hcop>
   Radio LTE Line Forest model <example_LTE>
   Near IR cube <example_NIR_cube>
   Optical SDSS <example_sdss>
   Optical Vega Echelle <example_vega_echelle>
   Optical INTERACTIVE <interactive>
   Optical Supernova (?) H-alpha multi-component voigt <example_sn>
   Optical MUSE Spectral Cube <example_MUSE>
   Template Fitting <example_template>
   Monte Carlo <example_pymc>
   Monte Carlo for NH3 <example_pymc_ammonia>

.. toctree::
   :maxdepth: 3
   :hidden:

   CLASS Guide <guide_class>
   IRAF Guide <guide_iraf>
   Student Projects <projects>

..
    TO DO:
    - add pymc / emcee examples
    - add UV examples?
    - add ALMA / EVLA cube examples


Additionally, here are some real-world examples.  They are not guaranteed to
work, but they show the use of pyspeckit in pipelines and paper-producing
code.

- `Fit each spectrum in a VLA NH3 6-6 cube (Goddi, Ginsburg, Zhang 2016) <https://github.com/keflavich/w51_vla_nh3_goddi/blob/master/code/fit_66_cube.py>`__
- `Fit each spectrum in an Arecibo H2CO 1-1 cube (Ginsburg et al, 2015) <https://github.com/keflavich/w51_singledish_h2co_maps/blob/master/analysis_scripts/pyspeckit_cube_fit_justlineprops.py>`__
- `Crop an ALMA CH3CN cube, then fit each spectrum <https://github.com/keflavich/SgrB2_ALMA_3mm_Mosaic/blob/9c3d4b79aa4ea01e848bd5c29c1949e9ba3b42e7/analysis/ch3cn_fiteach.py>`__
- `Fit each spectrum in an APEX H2CO 3-2 cube with 1 and 2 components (Ginsburg et al, 2016) <https://github.com/keflavich/APEX_CMZ_H2CO/blob/e298704dc395c3829a4b22118521d56912dd7c19/analysis/fit_the_brick.py#L31>`__
- `Compute moments, then fit each spectrum in a MUSE cube (McCleod et al, 2015) <https://github.com/keflavich/OrionNotebooks/blob/f890e10eaae479a95c766877d752bc38aca096f5/muse/vmap_085.py>`__
- `Fit H2 and BrG in a TripleSpec NIR cube (Youngblood et al, 2016) <https://github.com/allisony/TspecCubes/blob/1d02dbb8b274b0cc7fac0fcd378e2f2712434ff5/fit_cube.py>`__
- `Fit NH3 1-1 and 2-2 spectra for the GBT Ammonia Survey (Pineda, Friesen et al) <https://github.com/GBTAmmoniaSurvey/GAS/blob/9769e2e845915a2bb1801b66405a0f8eb787f5bc/GAS/PropertyMaps.py>`__
