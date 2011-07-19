import atpy
import urllib,urllib2
import tempfile

"""
Module to search Splatalogue.net via splat, modeled loosely on
ftp://ftp.cv.nrao.edu/NRAO-staff/bkent/slap/idl/
"""

length_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'angstroms':1e-10,'A':1e-10,
        }

# example query of SPLATALOGUE directly:
# http://www.cv.nrao.edu/php/splat/c.php?sid%5B%5D=64&sid%5B%5D=108&calcIn=&data_version=v2.0&from=&to=&frequency_units=MHz&energy_range_from=&energy_range_to=&lill=on&tran=&submit=Search&no_atmospheric=no_atmospheric&no_potential=no_potential&no_probable=no_probable&include_only_nrao=include_only_nrao&displayLovas=displayLovas&displaySLAIM=displaySLAIM&displayJPL=displayJPL&displayCDMS=displayCDMS&displayToyaMA=displayToyaMA&displayOSU=displayOSU&displayRecomb=displayRecomb&displayLisa=displayLisa&displayRFI=displayRFI&ls1=ls1&ls5=ls5&el1=el1

def query_splatalogue(minwav=0.00260,maxwav=0.00261,
        waveunits='m',root_url='http://find.nrao.edu/splata-slap/slap',
        chemical_element=None, NRAO_Recommended=True,):
    """
    Acquire an atpy table of a splatalogue searched based on wavelength.

    Future work will allow queries based on other parameters.  I'm waiting on
    development by the SPLATASLAP database folks to implement these.
    """

    if waveunits in length_dict:
        minwav = minwav * length_dict[waveunits]
        maxwav = maxwav * length_dict[waveunits]

    #query_url = "%s?REQUEST=queryData&WAVELENGTH=%f/%f" % (root_url,minwav,maxwav)
    # This is probably the more robust/pythonic way to do this sort of query:
    request_dict = {"REQUEST":"queryData","WAVELENGTH":"%f/%f" % (minwav,maxwav)}
    if chemical_element is not None: request_dict['CHEMICAL_ELEMENT'] = chemical_element
    #if NRAO_Recommended: request_dict['include_only_nrao'] = "yes"
    if NRAO_Recommended: request_dict['recommended'] = "1"
    #query_url = urllib2.Request(url=root_url,
    #        data=urllib.urlencode(request_dict))
    query_url = "%s?%s" % (root_url,urllib.urlencode(request_dict))

    U = urllib2.urlopen(query_url)

    # Normally it would be possible to do this:
    # t = atpy.Table(U,type='vo')
    # instead we have to write to a file and flush it. 
    # (see the error message below)

    R = U.read()
    U.close()
    tf = tempfile.NamedTemporaryFile()
    #for line in R:
    #    print >>tf,line.strip()
    print >>tf,R
    print tf.name
    tf.file.flush()
    t = atpy.Table(tf.name,type='vo')

    return t

"""
U = urllib2.urlopen(query_url)
t = atpy.Table(U,type='vo')
None:9:0: W35: 'value' attribute required for 'INFO' elements
None:18:0: W03: Implicitly generating an ID from a name 'catalog name' -> 'catalog_name'
None:27:0: W03: Implicitly generating an ID from a name 'molecular formula' -> 'molecular_formula'
None:33:0: W03: Implicitly generating an ID from a name 'molecule type' -> 'molecule_type'
None:54:0: W03: Implicitly generating an ID from a name 'frequency recommended' -> 'frequency_recommended'
None:57:0: W03: Implicitly generating an ID from a name 'quantum numbers' -> 'quantum_numbers'
------------------------------------------------------------
Traceback (most recent call last):
  File "<ipython console>", line 1, in <module>
  File "/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/atpy/basetable.py", line 167, in __init__
    self.read(*args, **kwargs)
  File "/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/atpy/basetable.py", line 213, in read
    atpy._readers[table_type](self, *args, **kwargs)
  File "/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/atpy/votable.py", line 86, in read
    votable = parse(filename, pedantic=pedantic)
  File "/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/vo/table.py", line 103, in parse
    return tree.VOTableFile(config=config, pos=(1, 1)).parse(iterator, config)
  File "/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/vo/tree.py", line 2921, in parse
    for start, tag, data, pos in iterator:
ValueError: 1:0: no element found
"""
