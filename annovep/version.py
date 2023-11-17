def _get_version() -> str:
    try:
        from annovep._version import VERSION
    except ModuleNotFoundError:
        return "Undefined"
    else:
        return VERSION


VERSION = _get_version()

__all__ = [
    "VERSION",
]
