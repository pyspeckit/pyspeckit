try:
    from astropy.models import ParametricModel,Parameter,_convert_input,_convert_output
    import numpy as np

    class PowerLawModel(ParametricModel):
        param_names = ['scale', 'alpha']

        def __init__(self, scale, alpha, param_dim=1):
            self._scale = Parameter(name='scale', val=scale, mclass=self, param_dim=param_dim)
            self._alpha = Parameter(name='alpha', val=alpha, mclass=self, param_dim=param_dim)
            super(ParametricModel,self).__init__(self, self.param_names, ndim=1, outdim=1, param_dim=param_dim)
            self.linear = False
            self.deriv = None
        
        def eval(self, xvals, params):
            return params[0]*((xvals)**(-params[1]))

        def noderiv(self, params, xvals, yvals):
            deriv_dict = {
                    'scale': ((xvals)**(-params[1])),
                    'alpha': params[0]*((xvals)**(-params[1]))*np.log(xvals)}
            derivval = [deriv_dict[par] for par in self.param_names]
            return np.array(derivval).T
            
        def __call__(self, x):
            """
            Transforms data using this model.
            
            Parameters
            --------------
            x : array, of minimum dimensions 1
            
            Notes
            -----
            See the module docstring for rules for model evaluation. 
            """
            x, fmt = _convert_input(x, self.param_dim)
            result = self.eval(x, self.param_sets)
            return _convert_output(result, fmt)

except ImportError:
    pass

