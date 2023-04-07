try:
    # Try getting the version without invoking setuptools-scm
    from ._version import version as __version__
except (ImportError, ModuleNotFoundError):
    try:
        from setuptools_scm import get_version
        __version__ = get_version()
    except ImportError:
        __version__ = "unknown"