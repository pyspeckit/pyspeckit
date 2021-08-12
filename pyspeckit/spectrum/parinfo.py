from __future__ import print_function
from six.moves import xrange
try:
    from lmfit import Parameters, Parameter
    LMFIT_PARAMETERS_INSTALLED = True
except ImportError:
    LMFIT_PARAMETERS_INSTALLED = False

class ParinfoList(list):
    """
    Store a list of model parameter values and their associated metadata (name,
    error, order, limits, etc.) in a class-friendly manner
    """
    def __init__(self, *args, **kwargs):
        """
        Create a parinfolist from either a list of parameters
        or a list or a tuple

        Parameters
        ----------
        preserve_order : bool
            Keep the p['n'] values for each parameter?
            Default to false, so par #'s will go from 0 to n(pars)-1
        """

        if LMFIT_PARAMETERS_INSTALLED:
            list.__init__(self,[])
            if len(args) == 1 and isinstance(args[0],Parameters):
                self._from_Parameters(args[0])
                self._dict = dict([(pp['parname'],pp) for pp in self])
                return

        list.__init__(self, *args)

        preserve_order = kwargs.pop('preserve_order',False)
        # re-order the parameters from 0 to n-1 unless told otherwise
        if not preserve_order:
            self._set_numbers()

        self._check_names()
        self._set_attributes()
        self._dict = dict([(pp['parname'],pp) for pp in self])

    def _set_numbers(self):
        """ Set the parameters in order by their current order in the list """
        for ii,pp in enumerate(self):
            if pp.n != ii:
                pp.n = ii

    def __str__(self):
        return "\n".join([repr(p) for p in self])

    def _getter(attributename):
        def getattribute(self):
            return [v[attributename] for v in self]
        return getattribute

    def _setter(attributename):
        def setattribute(self, values):
            if len(values) == len(self):
                for parinf,newval in zip(self,values):
                    parinf[attributename] = newval
            else:
                raise ValueError("Must have len(new values) = %i (was %i)" % (len(self),len(values)))
        return setattribute

    def keys(self):
        """ Dictionary-like behavior """
        return self.parnames

    def items(self):
        """ Dictionary-like behavior """
        return zip(self.parnames, self[:])

    #def values(self):
    #    """ Dictionary-like behavior """
    #    return [v['value'] for v in self]


    names = property(fget=_getter('parname'), fset=_setter('parname'))
    parnames=names
    shortnames = property(fget=_getter('shortparname'), fset=_setter('shortparname'))
    shortparnames=shortnames
    values = property(fget=_getter('value'), fset=_setter('value'))
    errors = property(fget=_getter('error'), fset=_setter('error'))
    n = property(fget=_getter('n'), fset=_setter('n'))
    order=n
    fixed = property(fget=_getter('fixed'), fset=_setter('fixed'))
    limits = property(fget=_getter('limits'), fset=_setter('limits'))
    limited = property(fget=_getter('limited'), fset=_setter('limited'))
    tied = property(fget=_getter('tied'), fset=_setter('tied'))

    def __getitem__(self, key):
        if type(key) in (int, slice):
            return super(ParinfoList,self).__getitem__(key)
        else:
            return self._dict[key]

    def __setitem__(self, key, val):
        """
        DO NOT allow items to be replaced/overwritten,
        instead use their own setters
        """
        # if key already exists, use its setter
        if key in self._dict or (type(key) is int and key < len(self)):
            self[key] = val
        elif type(key) is int:
            # can't set a new list element this way
            raise IndexError("Index %i out of range" % key)
        elif isinstance(val,Parinfo):
            # otherwise, add the item
            self.__dict__[key] = val
        else:
            raise TypeError("Can only add Parinfo items to ParinfoLists")

    def _set_attributes(self):
        self.__dict__.update(dict([(pp['parname'],pp) for pp in self]))

    def _check_names(self):
        """
        Make sure all names are unique.  If they're not, append #'s at the end
        (but strip #'s first)
        """
        name_counter = {}
        names_stripped = [name.strip('0123456789') for name in self.names]
        for ii,name in enumerate(names_stripped):
            if names_stripped.count(name) > 1:
                if name in name_counter:
                    name_counter[name] += 1
                    self[ii]['parname'] = self[ii]['parname'].strip('0123456789')+ "{0}".format(name_counter[name])
                else:
                    name_counter[name] = 0
                    self[ii]['parname'] = self[ii]['parname'].strip('0123456789')+ "{0}".format(name_counter[name])

                    # remove un-numbered versions if numbered versions are now being used
                    if name in self.__dict__:
                        self.__dict__.pop(name)

    def append(self, value, renumber=None):
        """
        Append to the list.  Will renumber the parameter if its number already
        exists in the list unless renumber == False
        """
        if hasattr(value,'n') and value.n in self.n and renumber is not False:
            # indexed from 0, so len(self) = max(self.n)+1
            value.n = len(self)
        super(ParinfoList, self).append(value)
        self._check_names()
        self._set_attributes()

    def as_Parameters(self):
        """
        Convert a ParinfoList to an lmfit Parameters class
        """
        if LMFIT_PARAMETERS_INSTALLED:
            P = Parameters()
            for par in self:
                P.add(name=par.parname,
                      value=par.value,
                      vary=not(par.fixed),
                      expr=par.tied if par.tied != '' else None,
                      min=par.limits[0] if par.limited[0] else None,
                      max=par.limits[1] if par.limited[1] else None)

            return P

    def _from_Parameters(self, lmpars):
        """
        Read from an lmfit Parameters instance
        """

        if len(lmpars) == len(self):
            for P in lmpars.values():
                self[P.name].value = P.value
                self[P.name].error = P.stderr
                self[P.name].limits = (P.min,P.max)
                self[P.name].limited = (P.min is not None,P.max is not None)
                self[P.name].expr = '' if P.expr is None else P.expr
        else:
            for par in lmpars.values():
                self.append(Parinfo(par))

    def tableprint(self, item_length=15, numbered=True):
        """
        Print data in table-friendly format

        Parameters
        ----------
        item_length : int
            Number of characters per item printed
        numbered : bool
            Are the parameters numbered?  In pyspeckit, they will always be,
            but this was included for compatibility with generic fitters
        """
        stripped_names = list(set([par.parname.strip("0123456789") for par in self]))

        nlines = len(self.n) / len(stripped_names)

        strformat = "%" + str(item_length) + "s"
        fltformat = "%" + str(item_length) + "g"

        print(" ".join([strformat % name for name in stripped_names]))
        if numbered:
            for ii in xrange(nlines):
                print(" ".join([fltformat % (self[name+"%i" % ii].value) for
                                name in stripped_names]))
        else:
            print(" ".join([fltformat % (self[name].value) for name in
                            stripped_names]))

class Parinfo(dict):
    """
    A class for storing attributes of a fitted model parameter.  It is based on
    mpfit's parinfo dictionary, which is just a dictionary containing a few set
    values.  This implements them as 'gettable' attributes instead, but far
    more importantly, includes sanity checks when setting values.

    Attributes
    ----------
    value: number
        The value of the parameter.  Arithmetic operations (*,/,+,-,**) will
        use this value
    error: number
        The error on the value
    n: int
        The order of the parameter in the model or ParinfoList
    fixed: bool
        Can the value change?  If False, error should be 0.
    limits: (min,max)
        The limits on the value of the parameter.  Only applied
        if limited
    limited: (bool, bool)
        Is the parameter value limited?
    step: number
        from MPFIT: the step size to be used in calculating the numerical
        derivatives.  If set to zero, then the step size is computed
        automatically.  Ignored when AUTODERIVATIVE=0.
    scaleable: bool
        Is the value scaled with the data?  Important for normalization
        procedures
    tied: string
        mpfit/lmift parameter.  Allows you to specify arbitrary expressions for
        how the parameter depends on other parameters
        mpfit:
            a string expression which "ties" the parameter to other free or
            fixed parameters.  Any expression involving constants and the
            parameter array P are permitted.  Example: if parameter 2 is always
            to be twice parameter 1 then use the following: parinfo(2).tied =
            '2 * p(1)'.  Since they are totally constrained, tied parameters
            are considered to be fixed; no errors are computed for them.
            NOTE: the PARNAME can't be used in expressions.
    parname: string
        The parameter name
    shortparname: string (tex)
        A shortened version of the parameter name for plotting display

    """
    def __init__(self, values=None, **kwargs):

        dict.__init__(self, {'value':0.0, 'error':0.0, 'n':0, 'fixed':False,
                             'limits':(0.0,0.0), 'limited':(False,False),
                             'step':False, 'scaleable':False, 'tied':'',
                             'parname':'', 'shortparname':''}, **kwargs)

        if LMFIT_PARAMETERS_INSTALLED:
            if isinstance(values,Parameter):
                self._from_Parameter(values)
                self.__dict__ = self
                return

        if values is not None:
            self.update(values)

        self.__dict__ = self

    def __repr__(self):
        try:
            reprint = "Param #%i %12s = %12g" % (self.n, self.parname, self.value)
            if self.fixed:
                reprint += " (fixed)"
            elif self.error is not None:
                reprint += " +/- %15g " % (self.error)
            if any(self.limited):
                lolim = "[%g," % self.limits[0] if self.limited[0] else "(-inf,"
                uplim = "%g]" % self.limits[1] if self.limited[1] else "inf)"
                myrange = lolim + uplim
                reprint += "  Range:%10s" % myrange
            if self.tied != '':
                reprint += " Tied: %s" % self.tied
            if self.shortparname != '':
                reprint += " Shortparname: %s" % self.shortparname

            return reprint
        except AttributeError:
            return super(Parinfo,self).__repr__()

    def __deepcopy__(self, memo):
        copy = Parinfo(self)
        copy.__dict__ = copy
        return copy

    def __copy__(self):
        copy = Parinfo(self)
        copy.__dict__ = copy
        return copy

    @property
    def max(self):
        return self.limits[1]

    @max.setter
    def max(self, value):
        self.limits = (self.limits[0], value)

    @property
    def min(self):
        return self.limits[0]

    @min.setter
    def min(self, value):
        self.limits = (value, self.limits[1])

    @property
    def vary(self):
        return not self.fixed

    @vary.setter
    def vary(self, value):
        self.fixed = not value

    @property
    def expr(self):
        return self.tied

    @expr.setter
    def expr(self, value):
        self._check_OK('tied',value)
        self.tied = value

    def __setattr__(self, key, value):
        # DEBUG print "Setting attribute %s = %s" % (key,value)
        self._check_OK(key,value)
        return super(Parinfo, self).__setattr__(key, value)

    def __setitem__(self, key, value):
        # DEBUG print "Setting item %s = %s" % (key,value)
        self._check_OK(key,value)
        return super(Parinfo, self).__setitem__(key, value)

    def _check_OK(self,key,value):
        # DEBUG print "Checking whether %s's value %s is OK" % (key,value)
        if key == 'value':
            if hasattr(self,'limited') and hasattr(self,'limits'):
                if self.limited[0] and value < self.limits[0]:
                    raise ValueError('Set parameter value %r < limit value %r' % (value,self.limits[0]))
                if self.limited[1] and value > self.limits[1]:
                    raise ValueError('Set parameter value %r > limit value %r' % (value,self.limits[1]))

        if key in ('limits','limited'):
            try:
                if len(value) != 2:
                    raise ValueError("%s must be a 2-tuple" % key)
            except TypeError: # if the input was scalar
                raise ValueError("%s must be a 2-tuple" % key)

        if key in ('parname','tied','shortparname'):
            if type(value) is not str:
                raise TypeError("%s must be a string" % key)

        if key in ('fixed',):
            try:
                value = bool(value)
            except:
                raise ValueError("%s had value %s, which could not be converted to boolean" % (key,value))

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

    def _from_Parameter(self, lmpar):
        """
        Read a Parinfo instance from and lmfit Parameter
        """
        self['limits'] = lmpar.min,lmpar.max
        self['limited'] = (lmpar.min not in (None,False),lmpar.max not in (None,False))
        self['value'] = lmpar.value
        self['error'] = lmpar.stderr
        self['parname'] = lmpar.name
        self['fixed'] = not(lmpar.vary)

    def _operation_wrapper(operation, reverse=False):
        """
        Perform an operation (addition, subtraction, mutiplication, division, etc.)
        """

        def ofunc(self, other):
            """ Operation Function """
            intypes = type(other),type(self.value)
            try:
                returnval = getattr(self.value,'__%s__' % operation)(other)
                if type(returnval) not in intypes:
                    raise TypeError("return value had wrong type: %s" % str(type(returnval)))
                else:
                    return returnval
            except TypeError as err: # integers don't have defined operations with floats
                #print err
                #print "TypeError1: ",self, self.value, other
                try:
                    if hasattr(other,'__r%s__' % operation):
                        #print "r",operation,": ",self, self.value, other
                        return getattr(other,'__r%s__' % operation)(self.value)
                    elif hasattr(other,'__%s__' % operation[1:]):
                        #print operation,": ",self, self.value, other
                        return getattr(other,'__%s__' % operation[1:])(self.value)
                except:
                    raise TypeError("Neither side of the operation has a %s attribute!" % operation)

        return ofunc

    __add__ = _operation_wrapper('add')
    __radd__ = _operation_wrapper('radd')
    __sub__ = _operation_wrapper('sub')
    __rsub__ = _operation_wrapper('rsub')
    __mul__ = _operation_wrapper('mul')
    __rmul__ = _operation_wrapper('rmul')
    __div__ = _operation_wrapper('div')
    __rdiv__ = _operation_wrapper('rdiv')
    __pow__ = _operation_wrapper('pow')
    __rpow__ = _operation_wrapper('rpow')

    try:
        def __array__(self):
            import numpy as np
            return np.array(self.value)
    except ImportError:
        pass


if __name__=="__main__":

    import unittest

    def check_failure(value):
        """ Raise a ValueError if value not in (0,1) """
        P = Parinfo()
        P.value = 1
        P.limited = (True,True)
        P.limits = (0,1)
        P.value = value

    def check_tied(value):
        P = Parinfo()
        P.tied = value

    def check_limits(value):
        P = Parinfo()
        P.limits = value

    def check_index(key):
        PL = ParinfoList([Parinfo({'parname':'HEIGHT'}),
                        Parinfo({'value':15,'parname':'AMPLITUDE'}),
                        Parinfo({'value':3,'parname':'WIDTH'}),
                        Parinfo({'value':4,'parname':'WIDTH'})])
        return PL[key]

    def check_set_list(values):
        PL = ParinfoList([Parinfo({'parname':'HEIGHT'}),
                        Parinfo({'value':15,'parname':'AMPLITUDE'}),
                        Parinfo({'value':3,'parname':'WIDTH','limits':(0,5),'limited':(True,True)}),
                        Parinfo({'value':4,'parname':'WIDTH'})])
        PL.shortparnames = ['a','b','c','d']
        PL.values = values
        return PL.values

    class MyTestCase(unittest.TestCase):
        def __init__(self, methodName='runTest', param=None):
            super(MyTestCase, self).__init__(methodName)
            self.param = param

        def test_checks_value_fail(self):
            check_failure(0.5)
            self.assertRaises(ValueError, check_failure, 5)
            self.assertRaises(ValueError, check_failure, -5)

        def test_checks_tied_fail(self):
            check_tied('p[0]')
            self.assertRaises(TypeError, check_tied, 5)
            self.assertRaises(TypeError, check_tied, (1,2,3))

        def test_checks_limits_fail(self):
            check_limits((1,2))
            self.assertRaises(ValueError, check_limits, -5)
            self.assertRaises(ValueError, check_limits, (1,2,3))

        def test_indexing(self):
            self.assertEqual(check_index(0), check_index('HEIGHT'))
            self.assertEqual(check_index(1), check_index('AMPLITUDE'))
            self.assertEqual(check_index(2), check_index('WIDTH0'))

        def test_set_list(self):
            self.assertEqual(check_set_list([1,2,3,4]),[1,2,3,4])
            self.assertRaises(ValueError,check_set_list,[1,2,10,4])

        def test_arithmetic(self):
            value = 25
            par = Parinfo({'parname':'TEST', 'value': value})
            assert par+5 == value+5
            assert par-5 == value-5
            assert par/5 == value/5
            assert par*5 == value*5

        def test_arithmetic2(self):
            value = 25.
            par = Parinfo({'parname':'TEST', 'value': value})
            assert par+5 == value+5
            assert par-5 == value-5
            assert par/5 == value/5
            assert par*5 == value*5

        def test_arithmetic3(self):
            value = 25.
            par = Parinfo({'parname':'TEST', 'value': value})
            assert par+5. == value+5.
            assert par-5. == value-5.
            assert par/5. == value/5.
            assert par*5. == value*5.

        def test_arithmetic4(self):
            value = 25
            par = Parinfo({'parname':'TEST', 'value': value})
            assert par+5. == value+5.
            assert par-5. == value-5.
            assert par/5. == value/5.
            assert par*5. == value*5.

        def test_arithmetic5(self):
            value = 25.
            par = Parinfo({'parname':'TEST', 'value': value})
            assert 5.+par == 5.+value
            assert 5.-par == 5.-value
            assert 5./par == 5./value
            assert 5.*par == 5.*value

        def test_arithmetic6(self):
            value = 25
            par = Parinfo({'parname':'TEST', 'value': value})
            assert 5.+par == 5.+value
            assert 5.-par == 5.-value
            assert 5./par == 5./value
            assert 5.*par == 5.*value

        def test_arithmetic7(self):
            value = 25.
            par = Parinfo({'parname':'TEST', 'value': value})
            assert 5+par == 5+value
            assert 5-par == 5-value
            assert 5/par == 5/value
            assert 5*par == 5*value

        def test_arithmetic8(self):
            value = 25
            par = Parinfo({'parname':'TEST', 'value': value})
            assert 5+par == 5+value
            assert 5-par == 5-value
            assert 5/par == 5/value
            assert 5*par == 5*value

        def test_copy(self):
            import copy
            value = 25
            par = Parinfo({'parname':'TEST', 'value': value})
            parcopy = copy.copy(par)
            assert parcopy.value == value
            parcopy = copy.deepcopy(par)
            assert parcopy.value == value

    unittest.main()
