mpfit_messages={
-16:"""A parameter or function value has become infinite or an undefined
   number.  This is usually a consequence of numerical overflow in the
   user's model function, which must be avoided.""",
0: "Improper input parameters.",
1: """Both actual and predicted relative reductions in the sum of squares
are at most ftol.""",
2: "Relative error between two consecutive iterates is at most xtol",
3: """Conditions for status = 1 and status = 2 both hold.""",
4: """The cosine of the angle between fvec and any column of the jacobian
  is at most gtol in absolute value.""",
5: """The maximum number of iterations has been reached.""",
6: """ftol is too small. No further reduction in the sum of squares is
  possible.""",
7: """xtol is too small. No further improvement in the approximate solution
  x is possible.""",
8: """gtol is too small. fvec is orthogonal to the columns of the jacobian
  to machine precision."""}

for ii in xrange(-15,0):
    mpfit_messages[ii] = """These are error codes that either MYFUNCT or iterfunct may return to
       terminate the fitting process.  Values from -15 to -1 are reserved
       for the user functions and will not clash with MPFIT."""

from ammonia import ammonia_model,ammonia_model_vtau
# old way from gaussfitter import gaussian_fitter
from inherited_gaussfitter import gaussian_fitter
from inherited_lorentzian import lorentzian_fitter
from inherited_voigtfitter import voigt_fitter
import formaldehyde
from formaldehyde import formaldehyde_fitter,formaldehyde_vheight_fitter
import n2hp
from n2hp import n2hp_vtau_fitter,n2hp_vtau
import hcn
import hyperfine,fitter,model,redshiftedgroup
import hill5infall
import radex_modelgrid
