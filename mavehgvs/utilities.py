from typing import List, Tuple, Optional

from mavehgvs.variant import Variant
from mavehgvs.exceptions import MaveHgvsParseError

__all__ = ["parse_variant_strings"]


def parse_variant_strings(
    variants: List[str]
) -> Tuple[List[Optional[Variant]], List[Optional[str]]]:
    """Parse a list of MAVE-HGVS strings into Variant objects or error messages.

    Parameters
    ----------
    variants : List[str]
        List of MAVE-HGVS strings to parse.

    Returns
    -------
    Tuple[List[Optional[Variant]], List[Optional[str]]]
        Returns a pair of lists containing variants or error messages.

        Both lists have the same length as the input list.
        The first list contains Variant objects if the string was successfully parsed; else None.
        The second list contains None if the string was successfully parsed; else the error message.

    """
    valid = list()
    invalid = list()

    for s in variants:
        try:
            v = Variant(s)
        except MaveHgvsParseError as error:
            valid.append(None)
            invalid.append(str(error))
        else:
            valid.append(v)
            invalid.append(None)

    return valid, invalid
