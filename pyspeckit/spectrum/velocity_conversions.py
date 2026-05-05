"""
Velocity-frame conversions (stub).

The original implementation depended on ``pyslalib``, which is not part of
the standard scientific Python stack any more.  Importing this module is a
no-op; calling :func:`convert_velocity` raises if the optional dependency
is missing.  Migrating these conversions to ``astropy.coordinates`` is
TODO.
"""
try:
    import pyslalib  # noqa: F401
    _PYSLALIB_AVAILABLE = True
except ImportError:
    _PYSLALIB_AVAILABLE = False


def convert_velocity(*args, **kwargs):
    """Placeholder for velocity-frame conversions; requires ``pyslalib``."""
    if not _PYSLALIB_AVAILABLE:
        raise ImportError("pyslalib is required for velocity frame conversions!")
    raise NotImplementedError(
        "velocity_conversions.convert_velocity is not implemented yet."
    )
