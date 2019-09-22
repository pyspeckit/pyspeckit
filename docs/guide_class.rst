Guide for GILDAS-CLASS users
============================

Reading Files
~~~~~~~~~~~~~
PySpecKit can read many file types including CLASS .cls files.

For example, a typical single-dish observing file will consist of a large number of 
individual spectra.  These may be on-the-fly spectra or multiple integrations on a single source;
either way the pointing and spectral information for each individual data point should be
incorporated in the data file.

The spectra are recorded from the data file into an
:class:`pyspeckit.classes.ObsBlock` object.   An ObsBlock (Observation Block) is simply a container
for many spectra; in most regards it behaves the same way as a :class:`pyspeckit.Spectrum` object,
but it has a few extra attributes built-in (e.g., spectral averaging).

A file can be loaded as a named object.
The data is filtered as it is read in.  The following line will ensure that the spectra in the
``n2hp`` ObsBlock contain only data from the F1M spectrometer with lines labeled 'N2HP(3-2)'
or 'N2H+(3-2)'::
    
    from pyspeckit.spectrum.readers.read_class import class_to_obsblocks
    n2hp = class_to_obsblocks(filename, 
        telescope=['SMT-F1M-HU','SMT-F1M-VU'],
        line=['N2HP(3-2)','N2H+(3-2)'])

Working with the data
~~~~~~~~~~~~~~~~~~~~~
You can treat an ObsBlock like any other `pyspeckit.Spectrum` object, but there
are a few unique features.

The 'smooth' function will smooth each individual spectrum in the ObsBlock, which 
can be useful if you want to smooth before averaging.

The 'average' function averages spectra.  Weights for the averaging can be specified.
The error is also computed either by taking the RMS of the averaged spectra or by
averaging the error spectra.

You can plot each spectrum in an individual window with the `ploteach` function, or all
overlaid simultaneously with the 'plotter' function.  If you want to loop through each
spectrum, waiting for user input after displaying, you can do something like::

    ax = pylab.gca()
    for sp in n2hp:
        sp.plotter(axis=ax)
        input("Waiting for input...")

You can also fit a line profile to *each* spectrum in the observation block
using the ``fiteach`` command.

Selecting Spectra
~~~~~~~~~~~~~~~~~

There are 3 keywords that can be used to select spectra when reading in a file.  The
``line`` and ``telescope`` keywords are required in order to make an
observation block, otherwise the spectral axes will not be common to all
spectra::

    n2hp = class_to_obsblocks(filename, 
        telescope=['SMT-F1M-HU','SMT-F1M-VU'],
        line=['N2HP(3-2)','N2H+(3-2)'])

If you want to read in all of the data and don't care about the line or telescope, you can instead use::

    all_data = class_to_spectra(filename)

You can select data after they are read in by matching header keywords::

    g10 = pyspeckit.ObsBlock([sp for sp in all_data if sp.specname == 'g10'])


Fitting Data
~~~~~~~~~~~~
In CLASS, you would specify the data range with individual commands on the command line, e.g.::

    [this isn't quite right]
    xrange 5 25

In pyspeckit, you can specify the range in multiple ways, but the default is to
use the plotted window.  For example::

    sp.plotter(xmin=5, xmax=25)
    sp.baseline()
    sp.specfit()

will perform a baseline fit and spectrum fit over the range 5-25 km/s (the
units of xmin, xmax are the plotted units).  More intricate specifications
are possible::

    sp.plotter()
    # Fit a baseline over the region 5-25 km/s, excluding 7-10 km/s and 15-18 km/s
    sp.baseline.selectregion(xmin=5,xmax=25,exclude=[7,10,15,18])
    sp.baseline()
    sp.specfit()
