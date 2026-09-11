import itertools
from typing import Iterable, Iterator, Tuple

from hypothesis import strategies as st

from ..strategies import (  # noqa: F401
    plain_positions,
    intron_offset_positions,
    intron_only_positions,
    utr_intron_positions,
    utr_only_positions,
    sequences,
)

# plain_positions / intron_offset_positions / intron_only_positions /
# utr_intron_positions / utr_only_positions / sequences are re-exported here
# (from tests/strategies.py, shared with tests/test_position.py and
# tests/test_variant.py) so that test_dna.py and test_rna.py can keep
# importing them from this package.


def nucleotide_variant_strings(  # noqa: max-complexity: 12
    position_strategy: st.SearchStrategy,
    alphabet: str,
    allow_bare_equal: bool = True,
) -> st.SearchStrategy:
    """Build a strategy generating single-variant strings (without a prefix) for
    a DNA- or RNA-style grammar: equal, sub, del, dup, ins, and delins events built
    from the given position strategy and nucleotide alphabet.

    This mirrors the six variant-event patterns combined into ``dna_variant_*``
    and ``rna_variant``.

    Set ``allow_bare_equal=False`` to exclude the bare ``=`` form of the equality
    event, guaranteeing the position strategy is actually used (useful for
    negative tests checking that a position is rejected).
    """

    @st.composite
    def equal(draw):
        kinds = ["position", "range"]
        if allow_bare_equal:
            kinds.append("bare")
        kind = draw(st.sampled_from(kinds))
        if kind == "bare":
            return "="
        elif kind == "position":
            return f"{draw(position_strategy)}="
        else:
            return f"{draw(position_strategy)}_{draw(position_strategy)}="

    @st.composite
    def sub(draw):
        pos = draw(position_strategy)
        ref = draw(st.sampled_from(alphabet))
        new = draw(st.sampled_from(alphabet))
        return f"{pos}{ref}>{new}"

    def _del_or_dup(suffix: str):
        @st.composite
        def strategy(draw):
            if draw(st.booleans()):
                return f"{draw(position_strategy)}{suffix}"
            else:
                return f"{draw(position_strategy)}_{draw(position_strategy)}{suffix}"

        return strategy()

    @st.composite
    def ins(draw):
        start = draw(position_strategy)
        end = draw(position_strategy)
        seq = draw(sequences(alphabet))
        return f"{start}_{end}ins{seq}"

    @st.composite
    def delins(draw):
        seq = draw(sequences(alphabet))
        if draw(st.booleans()):
            return f"{draw(position_strategy)}delins{seq}"
        else:
            start = draw(position_strategy)
            end = draw(position_strategy)
            return f"{start}_{end}delins{seq}"

    return st.one_of(
        equal(), sub(), _del_or_dup("del"), _del_or_dup("dup"), ins(), delins()
    )


def build_multi_variants(
    valid_strings: Iterable[str],
    invalid_strings: Iterable[str],
    min_length: int = 2,
    max_length: int = 3,
) -> Tuple[Iterator, Iterator]:
    """Build iterators of valid and invalid multi-variant strings to test.

    Parameters
    ----------
    valid_strings : Iterable[str]
        Iterable containing all the valid single-variant strings.
    invalid_strings : Iterable[str]
        Iterable containing all the invalid single-variant strings.
    min_length : int
        Minimum length of multi-variants that will be generated.
    max_length : int
        Maximum length of multi-variants that will be generated.
        Note that increasing this value may massively increase test runtime.

    Returns
    -------
    Tuple[Iterator, Iterator]
        Returns iterators containing semicolon-separated multi-variant strings.

        The first iterator contains multi-variants from only valid_strings and the
        second iterator contains multi-variants that include at least one variant from
        invalid_strings.
    """
    # create an iterable of permutations for each length and store them in lists
    valid_multivariants = list()
    invalid_multivariants = list()

    for i in range(min_length, max_length + 1):
        valid_multivariants.append(
            ";".join(x) for x in itertools.permutations(valid_strings, i)
        )
        invalid_multivariants.append(
            ";".join(x)
            for x in itertools.permutations(
                itertools.chain(valid_strings, invalid_strings), i
            )
            if any(y in x for y in invalid_strings)
        )

    # combine the lists into single iterators and return
    return itertools.chain.from_iterable(
        valid_multivariants
    ), itertools.chain.from_iterable(invalid_multivariants)
