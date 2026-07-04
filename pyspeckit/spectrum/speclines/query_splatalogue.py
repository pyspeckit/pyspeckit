"""
Deprecated Splatalogue SLAP query module.

The original implementation depended on ``atpy``, Python 2's ``urllib2``,
and the long-defunct NRAO "splata-slap" service; none of these are
available on modern stacks.  Use `astroquery.splatalogue
<https://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html>`_
instead.
"""

__all__ = ['query_splatalogue']


def query_splatalogue(*args, **kwargs):
    """
    Removed.  Use `astroquery.splatalogue` instead.
    """
    raise ImportError(
        "pyspeckit's legacy SLAP-based query_splatalogue has been removed "
        "(it depended on the unmaintained atpy package and a defunct web "
        "service).  Use astroquery.splatalogue instead: "
        "https://astroquery.readthedocs.io/en/latest/splatalogue/splatalogue.html"
    )
