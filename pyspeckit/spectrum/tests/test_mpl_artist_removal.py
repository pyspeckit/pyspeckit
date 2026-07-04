"""
Regression tests for issue #397.

matplotlib >= 3.7 turned ``Axes.lines``, ``Axes.patches``,
``Axes.collections``, and ``Axes.texts`` into immutable ``ArtistList``
views: calling ``.remove(artist)`` or ``.append(artist)`` on them raises
``AttributeError: 'ArtistList' object has no attribute 'remove'``.

These tests exercise the plot-artist removal paths (clearing fits,
re-fitting, baseline replotting, map click-marks/circles, line-ID
overlays, slider widgets) headlessly with the Agg backend and check that
the artists really are removed from the axes.
"""
import warnings

import matplotlib
matplotlib.use('Agg')

import numpy as np

from .. import Spectrum


def make_gaussian_spectrum():
    dx = 0.1
    x = np.arange(-6, 6, dx)
    y = np.exp(-x**2 / 2.) + 0.5
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=y, error=np.full(x.size, 0.1))
    return sp


def test_specfit_clear_removes_lines():
    """Specfit.clear() must not raise and must remove component lines."""
    sp = make_gaussian_spectrum()
    sp.plotter()
    sp.specfit(fittype='gaussian', guesses=[1, 0, 1])
    sp.specfit.plot_fit(show_components=True)
    axis = sp.plotter.axis

    assert len(sp.specfit._plotted_components) > 0
    n_before = len(axis.lines)

    sp.specfit.clear()

    assert sp.specfit._plotted_components == []
    assert len(axis.lines) < n_before


def test_refit_after_plotted_fit():
    """
    Reporter's scenario in #397: fitting a second time clears the plotted
    artifacts of the first fit and used to raise AttributeError.
    """
    sp = make_gaussian_spectrum()
    sp.plotter()
    sp.specfit(fittype='gaussian', guesses=[1, 0, 1])
    sp.specfit(fittype='gaussian', guesses=[1, 0, 1])
    # (the continuum offset biases the amplitude; we only care that the
    # second fit completed without an AttributeError from artist removal)
    assert np.isfinite(sp.specfit.parinfo[0].value)
    assert abs(sp.specfit.parinfo[1].value) < 0.1


def test_clear_highlights():
    """Interactive.clear_highlights hits axis.lines removal."""
    sp = make_gaussian_spectrum()
    sp.plotter()
    sp.specfit.selectregion(xmin=-3, xmax=3, highlight=True)
    axis = sp.plotter.axis

    assert len(sp.specfit.button1plot) > 0
    n_before = len(axis.lines)

    sp.specfit.clear_highlights()

    assert sp.specfit.button1plot == []
    assert len(axis.lines) < n_before


def test_baseline_overplot_twice():
    """
    Baseline.plotbaseline with subtract=False removes the previously
    plotted baseline lines on the second call.
    """
    sp = make_gaussian_spectrum()
    sp.plotter()
    sp.baseline(order=0, subtract=False)
    old_plots = list(sp.baseline._plots)
    assert len(old_plots) > 0

    sp.baseline(order=0, subtract=False)

    for p in old_plots:
        assert p not in sp.plotter.axis.lines


def test_baseline_subtract_replots():
    """
    Baseline with subtract=True clears every line from the axis before
    replotting the subtracted spectrum.
    """
    sp = make_gaussian_spectrum()
    # errstyle='fill' creates a PolyCollection errorplot, exercising the
    # errorplot-clearing branch of plotbaseline as well
    sp.plotter(errstyle='fill')
    sp.baseline(order=0, subtract=True)

    # the replotted (baseline-subtracted) spectrum must be there
    assert len(sp.plotter.axis.lines) > 0


def test_mapplot_click_marks_and_circles(tmp_path):
    """MapPlotter click marks and circles: add and remove."""
    from ...cubes.tests.test_cubetools import make_test_cube
    from ...cubes import Cube

    fname = str(tmp_path / 'test_mapplot.fits')
    make_test_cube((30, 9, 9), outfile=fname)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cube = Cube(fname)
    cube.mapplot(useaplpy=False, colorbar=False)
    mp = cube.mapplot

    mp._add_click_mark(3, 3)
    mp._add_click_mark(4, 4)
    assert len(mp._click_marks) == 2
    n_lines = len(mp.axis.lines)
    assert n_lines >= 2

    mp._clear_click_marks()
    assert mp._click_marks == []
    assert len(mp.axis.lines) == n_lines - 2

    mp._add_circle(3, 3, 5, 5)
    assert len(mp._circles) == 1
    assert len(mp.axis.patches) == 1

    mp._remove_circle()
    assert mp._circles == []
    assert len(mp.axis.patches) == 0


def test_radio_lines_hide():
    """radio_lines.hide removes texts and vline collections."""
    from ..speclines.radio import radio_lines

    sp = make_gaussian_spectrum()
    sp.plotter()
    axis = sp.plotter.axis

    # __init__ needs a splatalogue table, so build the instance manually
    # and populate it with real artists as show() would
    rl = radio_lines.__new__(radio_lines)
    rl.Spectrum = sp
    rl._lines = [axis.vlines([0.0, 1.0], 0, 1)]
    rl._linenames = [axis.text(0.0, 0.5, 'testline', rotation='vertical')]

    n_texts = len(axis.texts)
    n_collections = len(axis.collections)

    rl.hide()

    assert rl._lines == []
    assert rl._linenames == []
    assert len(axis.texts) == n_texts - 1
    assert len(axis.collections) == n_collections - 1


def test_modifiable_slider_valmin_valmax():
    """ModifiableSlider.set_valmin/set_valmax replace the vline in place."""
    import matplotlib.pyplot as plt
    from ..widgets import ModifiableSlider

    fig = plt.figure()
    ax = fig.add_subplot(111)
    slider = ModifiableSlider(ax, 'test', 0, 10, valinit=5)
    assert slider.vline in slider.ax.lines
    n_lines = len(slider.ax.lines)

    # valinit (5) < new valmin (6): triggers removal & re-creation of vline
    slider.set_valmin(6)
    assert slider.vline in slider.ax.lines
    assert len(slider.ax.lines) == n_lines

    # valinit (now 8) > new valmax (7): triggers the set_valmax branch
    slider.set_valmax(7)
    assert slider.vline in slider.ax.lines
    assert len(slider.ax.lines) == n_lines

    plt.close(fig)
