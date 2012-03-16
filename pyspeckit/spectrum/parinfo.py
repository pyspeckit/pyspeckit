
class ParinfoList(list):
    """
    Store a list of model parameter values and their associated metadata (name,
    error, order, limits, etc.) in a class-friendly manner
    """
    def __init__(self, *args):

        try:
            from lmfit import Parameters,Parameter
            list.__init__(self,[])
            if len(args) == 1 and isinstance(args[0],Parameters):
                self._from_Parameters(args[0])
                self._dict = dict([(pp['parname'],pp) for pp in self])
                return
        except ImportError:
            pass

        list.__init__(self, *args)

        for ii,pp in enumerate(self):
            if pp.n != ii:
                pp.n = ii

        self._check_names()
        self._set_attributes()
        self._dict = dict([(pp['parname'],pp) for pp in self])

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
        if type(key) is int:
            return super(ParinfoList,self).__getitem__(key)
        else:
            return self._dict[key]

    def _set_attributes(self):
        self.__dict__.update(dict([(pp['parname'],pp) for pp in self]))

    def _check_names(self):
        """
        Make sure all names are unique.  If they're not, append #'s at the end
        (but strip #'s first)
        """
        name_counter = {}
        names_stripped = [name.strip('[0-9]') for name in self.names]
        for ii,name in enumerate(names_stripped):
            if names_stripped.count(name) > 1:
                if name in name_counter:
                    name_counter[name] += 1
                    self[ii]['parname'] = self[ii]['parname'].strip('[0-9]')+ "{0}".format(name_counter[name])
                else:
                    name_counter[name] = 0
                    self[ii]['parname'] = self[ii]['parname'].strip('[0-9]')+ "{0}".format(name_counter[name])

                    # remove un-numbered versions if numbered versions are now being used
                    if name in self.__dict__:
                        self.__dict__.pop(name)

    def append(self, value):
        super(ParinfoList, self).append(value)
        self._check_names()
        self._set_attributes()

    def as_Parameters(self):
        """
        Convert a ParinfoList to an lmfit Parameters class
        """
        try:
            from lmfit import Parameters,Parameter
        except ImportError:
            print "Cannot import lmfit."
            return

        P = Parameters()
        P.add_many(*[(par.parname, par.value, not(par.fixed),
            par.limits[0] if par.limited[0] else None,
            par.limits[1] if par.limited[1] else None,
            par.tied if par.tied is not '' else None)
            for par in self])

        return P

    def _from_Parameters(self, lmpars):
        """
        Read from an lmfit Parameters instance
        """

        if len(lmpars) == len(self):
            self.names = [p.name for p in lmpars.values()]
            self.values = [p.value for p in lmpars.values()]
            self.errors = [p.stderr for p in lmpars.values()]
            for ii,P in enumerate(lmpars.values()):
                self[ii].limits = (P.min,P.max)
                self[ii].limited = (P.min not in (None,False),P.max not in (None,False))
                self[ii].expr = '' if P.expr is None else P.expr
        else:
            for par in lmpars.values():
                self.append(Parinfo(par))




class Parinfo(dict):
    """
    A class for storing attributes of a fitted model parameter.  It is based on
    mpfit's parinfo dictionary, which is just a dictionary containing a few set
    values.  This implements them as 'gettable' attributes instead, but far
    more importantly, includes sanity checks when setting values.
    """
    def __init__(self, values=None):

        dict.__init__(self, {'value':0.0, 'error':0.0,
                'n':0, 'fixed':False,
                'limits':(0.0,0.0),
                'limited':(False,False),
                'step':False,
                'tied':'',
                'parname':'',
                'shortparname':''})

        try:
            from lmfit import Parameter
            if isinstance(values,Parameter):
                self._from_Parameter(values)
                self.__dict__ = self
                return
        except ImportError:
            pass

        if values is not None:
            self.update(values)
        
        self.__dict__ = self

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

    unittest.main()
