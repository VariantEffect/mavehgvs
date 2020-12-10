__all__ = ["MaveHGVSParseError"]


class MaveHGVSParseError(Exception):
    """
    Base exception to use when a HGVS string cannot be parsed.
    """
    pass
