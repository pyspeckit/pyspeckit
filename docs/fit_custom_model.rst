Fitting a user-defined model to a spectrum
==========================================

``damped_lya_profile`` is the function that generates the model and ``damped_lya_fitter`` creates the fitter class.

.. code-block:: python

    def damped_lya_profile(wave_to_fit,vs_n,am_n,fw_n,vs_b,am_b,fw_b,h1_col,h1_b,
                           h1_vel,d2h,resolution,single_component_flux=False,
                           return_components=False,
                           return_hyperfine_components=False):

        """
        Computes a damped Lyman-alpha profile (by calling the functions
        lya_intrinsic_profile_func and total_tau_profile_func) and convolves it 
        to the proper resolution.

        """
    
        lya_intrinsic_profile = lya_intrinsic_profile_func(wave_to_fit,vs_n,am_n,fw_n,
                                 vs_b,am_b,fw_b,single_component_flux=single_component_flux)
        total_tau_profile = total_tau_profile_func(wave_to_fit,h1_col,h1_b,h1_vel,d2h)

        lya_obs_high = lya_intrinsic_profile * total_tau_profile

        ## Convolving the data ##
        fw = lya_rest/resolution
        aa = make_kernel(grid=wave_to_fit,fwhm=fw)
        lyman_fit = np.convolve(lya_obs_high,aa,mode='same')

        return lyman_fit*1e14

    def damped_lya_fitter(multisingle='multi'):
  
      """
      Generator for Damped LyA fitter class
  
      """
      myclass =  pyspeckit.models.model.SpectralModel(damped_lya_profile,
              11,
              parnames=['vs_n','am_n','fw_n','vs_b','am_b','fw_b','h1_col',
                        'h1_b','h1_vel','d2h','resolution'], 
              parlimited=[(False,False),(True,False),(True,False),(False,False),
                          (True,False),(True,False),(True,True),(True,True),
                          (False,False),(True,True),(True,False)], 
              parlimits=[(14.4,14.6), (1e-16,0), (50.,0),(0,0), (1e-16,0),(50.,0),
                         (17.0,19.5), (5.,20.), (0,0),(1.0e-5,2.5e-5),(10000.,0)])
      myclass.__name__ = "damped_lya"
      
      return myclass
  
    ## Load the spectrum, register the fitter, and fit
    spec = pyspeckit.Spectrum(xarr=wave_to_fit, data=flux_to_fit*1e14, 
                              error=error_to_fit*1e14, doplot=False,
                              header=spec_header)
  
  
    spec.Registry.add_fitter('damped_lya',damped_lya_fitter(),11)
  
    spec.specfit(fittype='damped_lya',
                 guesses=[vs_n, am_n, fw_n, vs_b, am_b, fw_b, h1_col, h1_b,
                          h1_vel, d2h, resolution],
                 quiet=False, fixed=[False, False, False, False, False,
                                     False, False, False, False, True, True])
  
    ## Check out the results
    spec.specfit.parinfo.values
    spec.specfit.parinfo.errors
    reduced_chi2 = spec.specfit.chi2/spec.specfit.dof
  
