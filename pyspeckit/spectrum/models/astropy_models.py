from astropy.models import ParametricModel,Parameter,_convert_input,_convert_output
import numpy as np

class PowerLawModel(ParametricModel):
    param_names = ['scale', 'alpha', 'characteristic']

    def __init__(self, scale, alpha, characteristic, param_dim=1):
        self.linear = False
        self._scale = Parameter(name='scale', val=scale, mclass=self, param_dim=param_dim)
        self._alpha = Parameter(name='alpha', val=alpha, mclass=self, param_dim=param_dim)
        self._characteristic = Parameter(name='characteristic', val=characteristic, mclass=self, param_dim=param_dim)
        ParametricModel.__init__(self, self.param_names, ndim=1, outdim=1, param_dim=param_dim)
    
    def eval(self, xvals, params):
        return params[0]*((xvals/params[2])**params[1])

    def deriv(self, params, xvals, yvals):
        return
        deriv_dict = {
                'scale': ((xvals/params[2])**params[1]),
                'alpha': params[0]*((xvals/params[2])**(params[1]-1))*np.log(1./params[2]),
                'characteristic': -1*params[0]*((xvals/params[2])**(params[1])) * (params[1]/params[2])}
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

