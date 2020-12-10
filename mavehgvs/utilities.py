from typing import List, Tuple, Optional

from mavehgvs.variant import Variant
from mavehgvs.exceptions import MaveHGVSParseError

__all__ = ["parse_variants"]


def parse_variants(variants: List[str]) -> Tuple[List[Optional[Variant]], List[Optional[str]]]:
    valid, invalid = [], []

    for variant in variants:
        try:
            v = Variant(variant)
            valid.append(v)
            invalid.append(None)
        except MaveHGVSParseError as error:
            valid.append(None)
            invalid.append(str(error))

    return valid, invalid
