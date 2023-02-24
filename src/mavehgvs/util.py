from typing import List, Tuple, Optional, Iterable

from mavehgvs.variant import Variant
from mavehgvs.exceptions import MaveHgvsParseError

__all__ = ["parse_variant_strings"]


def parse_variant_strings(
    variants: Iterable[str],
    targetseq: Optional[str] = None,
    expected_prefix: Optional[str] = None,
) -> Tuple[List[Optional[Variant]], List[Optional[str]]]:
    """Parse a list of MAVE-HGVS strings into Variant objects or error messages.

    Parameters
    ----------
    variants : Iterable[str]
        Iterable of MAVE-HGVS strings to parse.

    targetseq : Optional[str]
        If provided, all variants will be validated for agreement with this sequence.
        See the documentation for :py:class:`Variant` for further details.

    expected_prefix : Optional[str]
        If provided, all variants will be expected to have the same single-letter prefix.
        Variants that do not have this prefix will be treated as invalid.

    Returns
    -------
    Tuple[List[Optional[Variant]], List[Optional[str]]]
        Returns a pair of lists containing variants or error messages.

        Both lists have the same length as the input list.
        The first list contains Variant objects if the string was successfully parsed; else None.
        The second list contains None if the string was successfully parsed; else the error message.

    """
    if expected_prefix is not None and expected_prefix not in list("cgmnopr"):
        raise ValueError("invalid expected prefix")

    valid = list()
    invalid = list()

    for s in variants:
        try:
            v = Variant(s, targetseq=targetseq)
        except MaveHgvsParseError as error:
            valid.append(None)
            invalid.append(str(error))
        else:
            if expected_prefix is not None and v.prefix != expected_prefix:
                valid.append(None)
                invalid.append("unexpected variant prefix")
            else:
                valid.append(v)
                invalid.append(None)

    return valid, invalid
