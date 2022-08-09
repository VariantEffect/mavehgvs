import itertools
from typing import Iterable, Iterator, Tuple


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
