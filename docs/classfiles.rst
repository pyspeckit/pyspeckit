Gildas CLASS files
------------------

Pyspeckit is capable of reading files from some versions of CLASS.  The CLASS
developers have stated that the GILDAS file format is private and will remain so,
and therefore there are no guarantees that the CLASS reader will work for your 
file.

Nonetheless, if you want to develop in python instead of SIC, the :mod:`~pyspeckit.spectrum.readers.read_class`
module is probably the best way to access CLASS data.

The `CLASS file specification
<http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node56.html>`_ is incomplete,
so much of the data reading is hacked together.  The code style is based off of
Tom Robitaille's `idlsave <http://astrofrog.github.io/idlsave/>`_ package.

An example usage.  Note that ``telescope`` and ``line`` are NOT optional keyword arguments,
they are just specified as such for clarity ::

    n2hp = class_to_obsblocks(fn1, telescope=['SMT-F1M-HU','SMT-F1M-VU'],
        line=['N2HP(3-2)','N2H+(3-2)'])

This will generate a :class:`~pyspeckit.spectrum.ObsBlock` from all data tagged with the 'telescope'
flags listed and lines matching either of those above.  The data selection is equivalent to a combination of ::

    find /telescope SMT-F1M-HU
    find /telescope SMT-F1M-VU
    find /line N2HP(3-2)
    find /line N2H+(3-2)

ALL of the data matching those criteria will be included in an ObsBlock.  They
will then be accessible through the ObsBlock's `speclist` attribute, or just by
indexing the ObsBlock directly.

An essentially undocumented API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyspeckit.spectrum.readers.read_class
    :members:
