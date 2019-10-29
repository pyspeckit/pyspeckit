FITS
----

A minimal header should look like this::

    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                  -32 / array data type
    NAXIS   =                    2 / number of array dimensions
    NAXIS1  =                  659
    NAXIS2  =                    2
    CRPIX1  =                  1.0
    CRVAL1  =   -4953.029632560421
    CDELT1  =    212.5358581542998
    CTYPE1  = 'VRAD-LSR'
    CUNIT1  = 'm/s     '
    BUNIT   = 'K       '
    RESTFRQ =          110.20137E9
    SPECSYS = 'LSRK    '
    END

A fits file with a header as above should be easily readable without any user effort::

    sp = pyspeckit.Spectrum('test.fits')

If you have multiple spectroscopic axes, e.g. ::

    CRPIX1A =                  1.0
    CRVAL1A =    110.2031747948101
    CTYPE1A = 'FREQ-LSR'
    CUNIT1A = 'GHz     '
    RESTFRQA=            110.20137

you can load that axis with the 'wcstype' keyword::

    sp = pyspeckit.Spectrum('test.fits',wcstype='A')
    

If you have a .fits file with a non-linear X-axis that is stored in
the .fits file as data (as opposed to being implicitly included in
a heaer), you can load it using a custom .fits reader.
An example implementation is given in the `tspec_reader
<http://pyspeckit.bitbucket.org/path/to/tspec_reader.py>`_.
It can be registered using :doc:`registration`: ::

    tspec_reader = check_reader(tspec_reader.tspec_reader)
    pyspeckit.register_reader('tspec',tspec_reader,'fits')


API
~~~

.. automodule:: pyspeckit.spectrum.readers.fits_reader
    :members:
