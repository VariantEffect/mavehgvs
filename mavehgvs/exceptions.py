__all__ = ["MaveHgvsParseError"]


class MaveHgvsParseError(Exception):
    """
    Base exception to use when a HGVS string cannot be parsed.
    """
    pass
