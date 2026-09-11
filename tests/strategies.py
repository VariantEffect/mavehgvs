"""Shared Hypothesis strategies for building position, amino acid, and sequence
strings that match the grammars in mavehgvs.patterns.position and
mavehgvs.patterns.protein.

These are reused by tests/test_position.py (VariantPosition) and
tests/test_patterns/*.py (the regex patterns themselves), which both need to
generate strings conforming to the same numeric and amino-acid grammars.
"""

from hypothesis import strategies as st
from fqfa.constants import AA_CODES

AMINO_ACIDS = sorted(set(AA_CODES.values()))


@st.composite
def plain_positions(draw) -> str:
    """A plain position with no intron or UTR component, e.g. '8'. Matches the
    ``pos`` pattern."""
    return str(draw(st.integers(min_value=1, max_value=10**9)))


@st.composite
def intron_offset_strings(draw) -> str:
    """An intronic offset with sign, e.g. '+10' or '-6'."""
    sign = draw(st.sampled_from(["+", "-"]))
    n = draw(plain_positions())
    return f"{sign}{n}"


@st.composite
def utr_position_strings(draw) -> str:
    """A UTR position with no intronic offset, e.g. '*8' (3') or '-80' (5')."""
    n = draw(plain_positions())
    if draw(st.booleans()):
        return f"*{n}"
    else:
        return f"-{n}"


@st.composite
def intron_offset_positions(draw) -> str:
    """A position with an optional intronic offset, e.g. '8' or '78+10'. Matches
    the ``pos_intron`` pattern."""
    base = draw(plain_positions())
    if draw(st.booleans()):
        base += draw(intron_offset_strings())
    return base


@st.composite
def intron_only_positions(draw) -> str:
    """A position guaranteed to include an intronic offset, e.g. '78+10'."""
    return draw(plain_positions()) + draw(intron_offset_strings())


@st.composite
def utr_intron_positions(draw) -> str:
    """A position with an optional UTR prefix and/or intronic offset, e.g. '*8',
    '-80', or '-45-1'. Matches the ``pos_intron_utr`` pattern."""
    base = draw(st.one_of(plain_positions(), utr_position_strings()))
    if draw(st.booleans()):
        base += draw(intron_offset_strings())
    return base


@st.composite
def utr_only_positions(draw) -> str:
    """A position guaranteed to include a UTR prefix, e.g. '*8' or '-45-1'."""
    base = draw(utr_position_strings())
    if draw(st.booleans()):
        base += draw(intron_offset_strings())
    return base


@st.composite
def amino_acid_position_strings(draw) -> str:
    """A protein position, e.g. 'Gly8'. Matches the ``aa_pos`` pattern. Amino
    acid positions cannot use the extended (UTR/intron) syntax."""
    aa = draw(st.sampled_from(AMINO_ACIDS))
    n = draw(plain_positions())
    return f"{aa}{n}"


@st.composite
def amino_acid_sequences(draw, min_size: int = 1, max_size: int = 5) -> str:
    """A sequence of concatenated three-letter amino acid codes, e.g. 'GlyProCys'.
    Used for protein insertions and deletion-insertions."""
    n = draw(st.integers(min_value=min_size, max_value=max_size))
    return "".join(draw(st.sampled_from(AMINO_ACIDS)) for _ in range(n))


def sequences(
    alphabet: str, min_size: int = 1, max_size: int = 100
) -> st.SearchStrategy:
    """A strategy for strings drawn from the given alphabet, e.g. nucleotide
    sequences for insertions and deletion-insertions."""
    return st.text(alphabet=alphabet, min_size=min_size, max_size=max_size)
