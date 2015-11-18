from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits

def intersection(header1, header2, if_conflict=None):
    """
    Return a pyfits Header containing the intersection of two pyfits Headers

    *if_conflict* [ '1'/1/'Header1' | '2'/2/'Header2' | None ]
        Defines behavior if a keyword conflict is found.  Default is to remove the key

    """

    newheader = pyfits.Header()

    for key,value in header1.items():
        if key in header2:
            try:
                if value == header2[key]:
                    newheader[key] = value
                elif if_conflict in ('1',1,'Header1'):
                    newheader[key] = value
                elif if_conflict in ('2',2,'Header2'):
                    newheader[key] = Header2[key]
            except KeyError:
                """ Assume pyfits doesn't want you to have that keyword
                (because it shouldn't be possible to get here otherwise) """
                pass
        else:
            try:
                newheader[key] = value
            except KeyError:
                """ Assume pyfits doesn't want you to have that keyword
                (because it shouldn't be possible to get here otherwise) """
                pass

    return newheader
