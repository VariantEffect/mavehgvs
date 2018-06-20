"""
HGVS regex parsing for DNA, RNA and Protein specifications. `Event` refers
to a specific mutation event which does not contain any prefixes i.e.
a character from 'pcngmr' prefixing a mutation syntax like p.Leu4Gly. A
`variant` refers to an `Event` prefixed with a the former characters, which
may also be enclosed in parentheses or brackets e.g. p.(Leu5Gly)
or p.[Leu4Gly;Gly7Leu]
"""

import re
from enum import Enum

from .dna import single_variant as dna_single_variant
from .dna import multi_variant as dna_multi_variant

from .rna import single_variant as rna_single_variant
from .rna import multi_variant as rna_multi_variant

from .protein import single_variant as protein_single_variant
from .protein import multi_variant as protein_multi_variant


class Event(Enum):
    """
    Enum for supported mutation events.
    """
    DELETION = 'del'
    INSERTION = 'ins'
    DELINS = 'delins'
    SUBSTITUTION = 'sub'
    FRAME_SHIFT = 'fs'
    
    @classmethod
    def str_to_enum(cls, item):
        if item[0] in ['>', '=', '=//', '=/']:
            return cls.SUBSTITUTION
        else:
            return cls._value2member_map_.get(item, None)


class Level(Enum):
    """
    Enum specifying supported HGVS syntax types.
    """
    DNA = 'dna'
    RNA = 'rna'
    PROTEIN = 'protein'
    
    @classmethod
    def str_to_enum(cls, item):
        return cls._value2member_map_.get(item, None)


# Remove capture groups used for use in joining regexes in
# multi-variants since capture groups cannot be defined more than once.
single_variant = r"({0})|({1})|({2})".format(
    dna_single_variant,
    rna_single_variant,
    protein_single_variant,
)
multi_variant = r"({0})|({1})|({2})".format(
    dna_multi_variant,
    rna_multi_variant,
    protein_multi_variant,
)


# ---- Compiled Regex Expressions
single_variant_re = re.compile(single_variant)
multi_variant_re = re.compile(multi_variant)


def infer_type(hgvs):
    """
    Infer the type of hgvs variant as supported by `Event`. The order of
    this if/elif block is important: DELINS should come before INSERTION and
    DELETION since DELINS events contain INSERTION and DELETION event values
    as substrings.
    
    Parameters
    ----------
    hgvs : str
        The hgvs string to infer from.
    
    Returns
    -------
    `Event`
        An `Enum` value from the `Event` enum.
    """
    if not hgvs:
        return None
    if Event.DELINS.value in hgvs:
        return Event.DELINS
    elif Event.INSERTION.value in hgvs:
        return Event.INSERTION
    elif Event.DELETION.value in hgvs:
        return Event.DELETION
    elif Event.FRAME_SHIFT.value in hgvs:
        return Event.FRAME_SHIFT
    else:
        return Event.SUBSTITUTION
    
    
def infer_level(hgvs):
    """
    Infer the level of hgvs variant as supported by `Level` by inspecting
    the prefix.

    Parameters
    ----------
    hgvs : str
        The hgvs string to infer from.

    Returns
    -------
    `Level`
        An `Enum` value from the `Level` enum.
    """
    if not hgvs:
        return None
    if hgvs[0] in 'cgnm':
        return Level.DNA
    elif hgvs[0] == 'r':
        return Level.RNA
    elif hgvs[0] == 'p':
        return Level.PROTEIN
    else:
        return None
    

def is_multi(hgvs):
    return bool(multi_variant_re.fullmatch(hgvs))
