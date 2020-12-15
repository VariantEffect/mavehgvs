from typing import List, Tuple, Optional

from mavehgvs.variant import Variant
from mavehgvs.exceptions import MaveHgvsParseError

__all__ = ["parse_variant_strings"]


def parse_variant_strings(
    variants: List[str]
) -> Tuple[List[Optional[Variant]], List[Optional[str]]]:
    valid, invalid = [], []

    for variant in variants:
        try:
            v = Variant(variant)
            valid.append(v)
            invalid.append(None)
        except MaveHgvsParseError as error:
            valid.append(None)
            invalid.append(str(error))

    return valid, invalid
