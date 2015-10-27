from . import classes
from .utils import reddening

class OpticalSpectrum(classes.Spectrum):
    def deredden(self,**kwargs):
        """ 
        Deredden a spectrum inplace

        Parameters
        ----------
        ebv: float
            E(B-V) differential extinction; specify either this or a_v.
        a_v: float
            A(V) extinction; specify either this or ebv.
        r_v: float, optional
            defaults to standard Milky Way average of 3.1
        model: {'f99', 'fm07'}, optional
            * 'f99' is the default Fitzpatrick (1999) [1]_
            * 'fm07' is Fitzpatrick & Massa (2007) [2]_. Currently not R dependent.
        model: {'ccm89', 'gcc09'}, optional
            * 'ccm89' is the default Cardelli, Clayton, & Mathis (1989) [1]_, but 
              does include the O'Donnell (1994) parameters to match IDL astrolib.
            * 'gcc09' is Gordon, Cartledge, & Clayton (2009) [2]_. This paper has
              incorrect parameters for the 2175A bump; not yet corrected here.
        """

        self.data *= self._get_reddening(**kwargs)

    def _get_reddening(self,ebv=None,a_v=None,r_v=3.1,model='ccm89'):

        if model in ('ccm89','gcc09'):
            dered = reddening.ccm_reddening(self.xarr.as_unit('angstrom'),
                    ebv=ebv, a_v=a_v, r_v=r_v, model=model)
        elif model in ('f99','fm07'):
            dered = reddening.fm_reddening(self.xarr.as_unit('angstrom'),
                    ebv=ebv, a_v=a_v, r_v=r_v, model=model)

        return dered


    def redden(self,**kwargs):
        """ Un-de-redden a spectrum in-place (see `deredden`) """

        self.data /= self._get_reddening(**kwargs)

