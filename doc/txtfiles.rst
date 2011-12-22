Plain Text
-------------------------------------

Text files should be of the form::

    wavelength flux err 
    3637.390 0.314 0.000
    3638.227 0.717 0.000
    3639.065 1.482 0.000

where there 'err' column is optional but the others are not.  The most basic
spectrum file allowed would have no header and two columns, e.g.::

    1  0.5
    2  1.5
    3  0.1

If the X-axis is not monotonic, the data will be sorted so that the X-axis is
in ascending order.


API
~~~

.. automodule:: pyspeckit.spectrum.readers.txt_reader
    :members:
