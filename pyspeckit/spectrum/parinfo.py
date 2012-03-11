
class ParinfoList(list):
    def __init__(self, *args):
        list.__init__(self, *args)

        for ii,pp in enumerate(self):
            if pp.n != ii:
                pp.n = ii

        self._check_names()
        self._set_attributes()

    @property
    def names(self):
        return [v.parname for v in self]

    parnames=names

    @property
    def values(self):
        return [v.value for v in self]

    @property
    def errors(self):
        return [v.error for v in self]


    @property
    def n(self):
        return [v.n for v in self]

    order=n

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



class Parinfo(dict):
    def __init__(self, values=None):
        dict.__init__(self, {'value':0.0, 'error':0.0,
                'n':0, 'fixed':(False,False),
                'limits':(0.0,0.0),
                'limited':(False,False),
                'step':False,
                'tied':'',
                'parname':''})

        if values is not None:
            self.update(values)
        
        self.__dict__ = self

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

        if key in ('limits','limited','fixed'):
            try:
                if len(value) != 2:
                    raise ValueError("%s must be a 2-tuple" % key)
            except TypeError: # if the input was scalar
                raise ValueError("%s must be a 2-tuple" % key)

        if key in ('parname','tied'):
            if type(value) is not str:
                raise TypeError("%s must be a string" % key)

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
    
    class MyTestCase(unittest.TestCase):
        def test_checks(self):
            check_failure(0.5)
            self.assertRaises(ValueError, check_failure, 5)
            self.assertRaises(ValueError, check_failure, -5)
            check_tied('p[0]')
            self.assertRaises(TypeError, check_tied, 5)
            self.assertRaises(TypeError, check_tied, (1,2,3))
            check_limits((1,2))
            self.assertRaises(ValueError, check_limits, -5)
            self.assertRaises(ValueError, check_limits, (1,2,3))

    unittest.main()
