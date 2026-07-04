"""
Regression tests for the atpy -> astropy.table migration (issue #325):
speclines.radio.get_splat_table and Spectra.fiteach's fittable.
"""
import warnings

import numpy as np
import pytest
from astropy.table import Table

from .. import Spectrum, Spectra
from ..speclines import radio

# Header (and format) of the historical colon-delimited splatalogue.csv
# export that get_splat_table was written to parse.  The first data row has
# an empty computed frequency to exercise the masked-value handling.
SPLAT_CSV = """Species:Chemical Name:Freq-GHz:Freq Err:Meas Freq-GHz:Meas Freq Err:Resolved QNs:CDMS/JPL Intensity:Lovas/AST<br>Intensity:E<sub>L</sub> (K):E<sub>U</sub> (K):Linelist
t-HCOOH:Formic Acid:::0.00035:3.0E-7:3(3,0)-3(3,1):-15.05800::35.09214:35.09216:CDMS
O2v=0:Molecular oxygen:0.00050:1.0E-7:::N=3,J=3:-17.85120::26.38299:26.38301:JPL
O2v=0:Molecular oxygen:0.00100:1.0E-7:::N=1,J=1:-16.81010::5.69911:5.69915:JPL
HC7Nv=0:Cyanohexatriyne:23.6879:0.0001:23.6879:0.0001:J=21-20:-3.46490::11.50465:12.64150:CDMS
"""


@pytest.fixture
def splat_csv(tmp_path):
    fn = tmp_path / "splatalogue.csv"
    fn.write_text(SPLAT_CSV)
    return str(fn)


def test_get_splat_table_astropy(splat_csv):
    splat = radio.get_splat_table(filename=splat_csv)

    assert isinstance(splat, Table)

    # column names are stripped of non-alphanumeric characters
    for colname in ('Species', 'ChemicalName', 'FreqGHz', 'MeasFreqGHz',
                    'ResolvedQNs', 'Linelist'):
        assert colname in splat.colnames

    # columns synthesized by get_splat_table
    for colname in ('LineName', 'LatexName', 'frequency'):
        assert colname in splat.colnames

    # frequency columns must be float, with missing values filled with -999
    assert splat['FreqGHz'].dtype.kind == 'f'
    assert splat['MeasFreqGHz'].dtype.kind == 'f'
    assert not np.any(np.isnan(splat['FreqGHz']))
    np.testing.assert_allclose(np.asarray(splat['FreqGHz']),
                               [-999, 0.0005, 0.001, 23.6879])
    np.testing.assert_allclose(np.asarray(splat['frequency']),
                               np.asarray(splat['FreqGHz']))
    # MeasFreqGHz missing in rows 2 & 3
    np.testing.assert_allclose(np.asarray(splat['MeasFreqGHz']),
                               [0.00035, -999, -999, 23.6879])

    # LineName/LatexName are built from Species + ResolvedQNs
    assert splat['LineName'][3] == 'HC7Nv=0J=21-20'


def test_get_splat_table_missing_file(tmp_path):
    with pytest.raises(IOError):
        radio.get_splat_table(filename=str(tmp_path / "does_not_exist.csv"))


def test_get_splat_table_webquery_deprecated(splat_csv):
    with pytest.warns(DeprecationWarning):
        splat = radio.get_splat_table(webquery=True, filename=splat_csv)
    assert len(splat) == 4


def test_spectra_fiteach_fittable():
    """Spectra.fiteach should build its fit-results table with astropy."""
    dx = 0.1
    x = np.arange(-6, 6, dx)
    y = np.exp(-x**2 / 2.)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sp1 = Spectrum(xarr=x, data=y, xarrkwargs={'unit': 'km/s'})
        sp2 = Spectrum(xarr=x, data=2 * y, xarrkwargs={'unit': 'km/s'})
        spectra = Spectra([sp1, sp2])

    spectra.fiteach(fittype='gaussian', guesses=[1.5, 0, 1])

    assert isinstance(spectra.fittable, Table)
    for colname in ('name', 'amplitude', 'center', 'width',
                    'amplitudeerr', 'centererr', 'widtherr'):
        assert colname in spectra.fittable.colnames
    assert len(spectra.fittable) == 2
    np.testing.assert_allclose(spectra.fittable['amplitude'], [1, 2],
                               rtol=1e-3)
