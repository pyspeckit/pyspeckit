Registration
------------

PySpecKit is made extensible by allowing user-registered modules for reading,
writing, and fitting data.

For examples of registration in use, look at the source code of
``pyspeckit.spectrum.__init__`` and  ``pyspeckit.spectrum.fitters``.


The registration functions can be accessed directly: ::

    pyspeckit.register_reader
    pyspeckit.register_writer

However, models are bound to individual instances of the Spectrum class, so they
must be accessed via a specfit instance ::

    sp = pyspeckit.Spectrum('myfile.fits')
    sp.specfit.register_fitter

Alternatively, you can access and edit the default Registry ::

    pyspeckit.fitters.default_Registry.add_fitter

If you've already loaded a Spectrum instance, but then you want to reload fitters
from the default_Registry, or if you want to make your own `Registry`, you can use
the semi-private method ::

    MyRegistry = pyspeckit.fitters.Registry()
    sp._register_fitters(registry=MyRegistry)


Examples
--------

If you want to register a new variable-optical-depth deuterated ammonia model,
you could do the following:

    sp.specfit.register_fitter(name='nh2d', function=nh2d.nh2d_vtau_fitter, npars=4)
    

API
~~~

.. automodule:: pyspeckit.spectrum.__init__
    :members:

.. automodule:: pyspeckit.spectrum.fitters
    :members:
