"""
Regression tests for issue #414: on matplotlib >= 3.8,
``FigureManagerBase.key_press`` no longer exists, so starting an interactive
fitting session printed

    Error '...' object has no attribute 'key_press' was raised when trying to
    connect the key_press handler.  Please file an issue on github. ...

and matplotlib's built-in key-press handler was never actually disconnected.
"""
import warnings

import numpy as np

from .. import Spectrum


def _make_spectrum():
    x = np.linspace(-50, 50, 101)
    y = np.exp(-x**2 / (2 * 5.0**2))
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        return Spectrum(xarr=x, data=y, xarrkwargs={'unit': 'km/s'})


def test_interactive_session_no_keypress_error(capsys):
    """
    Starting and ending an interactive fit session must not emit the
    "Please file an issue on github" key_press error, and must
    disconnect/restore matplotlib's default key handler.
    """
    sp = _make_spectrum()
    sp.plotter()

    manager = sp.plotter.figure.canvas.manager
    original_handler_id = getattr(manager, 'key_press_handler_id', None)

    sp.specfit(interactive=True)

    if original_handler_id is not None:
        # matplotlib's default key handler must be disconnected while the
        # interactive session is active, so that its keybindings ('g', 'l',
        # 'k', ...) do not conflict with pyspeckit's
        assert manager.key_press_handler_id is None

    # end the interactive session (equivalent to pressing 'd' = done)
    sp.specfit.clear_all_connections()

    if original_handler_id is not None:
        # the default handler must be restored afterwards
        assert manager.key_press_handler_id is not None

    out, err = capsys.readouterr()
    combined = out + err
    assert 'Please file an issue on github' not in combined
    assert 'was raised when trying to connect the key_press handler' not in combined

    # plotting and fitting must still work after the interactive session
    sp.plotter()
    sp.specfit(fittype='gaussian', guesses=[1, 0, 5])
    np.testing.assert_allclose(sp.specfit.parinfo.values, [1, 0, 5], atol=1e-3)


def test_disconnect_reconnect_idempotent():
    """
    The disconnect/reconnect helpers must be safe to call repeatedly and in
    either order (e.g. clear_all_connections calls reconnect before any
    disconnect has ever happened).
    """
    sp = _make_spectrum()
    sp.plotter()

    manager = sp.plotter.figure.canvas.manager
    original_handler_id = getattr(manager, 'key_press_handler_id', None)

    # reconnect with nothing disconnected: no-op, no error
    sp.plotter._reconnect_matplotlib_keys()
    assert getattr(manager, 'key_press_handler_id', None) == original_handler_id

    sp.plotter._disconnect_matplotlib_keys()
    sp.plotter._disconnect_matplotlib_keys()
    sp.plotter._reconnect_matplotlib_keys()
    sp.plotter._reconnect_matplotlib_keys()

    if original_handler_id is not None:
        assert manager.key_press_handler_id is not None
