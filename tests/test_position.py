import unittest
import itertools
import random
import string

from hypothesis import given, strategies as st

from mavehgvs.position import VariantPosition
from mavehgvs.exceptions import MaveHgvsParseError
from .strategies import (
    plain_positions,
    utr_position_strings,
    utr_intron_positions,
    amino_acid_position_strings,
)

# strategies for building position strings that are valid according to the
# grammar in mavehgvs.patterns.position / mavehgvs.patterns.protein; the
# building blocks live in tests/strategies.py, shared with tests/test_patterns/


@st.composite
def valid_position_strings(draw) -> str:
    """Any valid position: nucleotide or amino acid."""
    return draw(st.one_of(utr_intron_positions(), amino_acid_position_strings()))


# strategies for building pairs of position strings that are adjacent by
# construction, mirroring the categories of adjacent_pairs in
# TestAdjacency.test_adjacent_pairs


@st.composite
def sequential_same_kind_pairs(draw) -> tuple[str, str]:
    """Two non-intronic positions of the same kind whose numeric positions
    differ by exactly one, e.g. ('8', '9') or ('*1', '*2')."""
    kind = draw(st.sampled_from(["plain", "five_prime_utr", "three_prime_utr"]))
    n = draw(st.integers(min_value=1, max_value=10**9))
    if kind == "plain":
        return str(n), str(n + 1)
    elif kind == "five_prime_utr":
        return f"-{n + 1}", f"-{n}"
    else:
        return f"*{n}", f"*{n + 1}"


@st.composite
def adjacent_intronic_offset_pairs(draw) -> tuple[str, str]:
    """Two positions with the same base and intronic offsets of the same sign
    that differ by exactly one, e.g. ('99+88', '99+89') or ('100-12', '100-11')."""
    base = draw(st.one_of(plain_positions(), utr_position_strings()))
    sign = draw(st.sampled_from(["+", "-"]))
    n = draw(st.integers(min_value=1, max_value=10**9))
    return f"{base}{sign}{n}", f"{base}{sign}{n + 1}"


@st.composite
def intron_exon_boundary_pairs(draw) -> tuple[str, str]:
    """A base position paired with the first base of its adjacent intron, e.g.
    ('202-1', '202') or ('99', '99+1')."""
    base = draw(st.one_of(plain_positions(), utr_position_strings()))
    sign = draw(st.sampled_from(["+", "-"]))
    return base, f"{base}{sign}1"


@st.composite
def adjacent_position_pairs(draw) -> tuple[str, str]:
    """A pair of position strings that are adjacent by construction.

    Covers the same categories as the fixed examples in
    TestAdjacency.test_adjacent_pairs: sequential positions of the same kind,
    sequential intronic offsets, the exon/intron boundary, and the 5' UTR /
    coding sequence boundary at -1/1.
    """
    return draw(
        st.one_of(
            sequential_same_kind_pairs(),
            adjacent_intronic_offset_pairs(),
            intron_exon_boundary_pairs(),
            st.just(("-1", "1")),
        )
    )


class TestObjectCreation(unittest.TestCase):
    def test_position_only(self) -> None:
        v = VariantPosition("8")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (8, None, None, None),
        )
        self.assertFalse(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertFalse(v.is_extended())

        v = VariantPosition("92380")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (92380, None, None, None),
        )
        self.assertFalse(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertFalse(v.is_extended())

    def test_amino_acid(self) -> None:
        v = VariantPosition("Gly8")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (8, "Gly", None, None),
        )
        self.assertFalse(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertTrue(v.is_protein())
        self.assertFalse(v.is_extended())

        v = VariantPosition("Cys92380")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (92380, "Cys", None, None),
        )
        self.assertFalse(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertTrue(v.is_protein())
        self.assertFalse(v.is_extended())

    def test_invalid_strings(self) -> None:
        position_strings = (
            "08",
            "+12",
            "*-99",
            "A",
            "TCGA",
            "g",
            "*",
            "-",
            "+",
            "**6",
            "800 + 12",
            "-12*5",
            "Glu-12",
            "*5Trp",
            "Xyz12",
            "ALA12",
        )
        for s in position_strings:
            with self.subTest(s=s):
                with self.assertRaises(MaveHgvsParseError):
                    VariantPosition(s)

    def test_utr(self) -> None:
        v = VariantPosition("*8")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (8, None, None, True),
        )
        self.assertTrue(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

        v = VariantPosition("-80")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-80, None, None, True),
        )
        self.assertTrue(v.is_utr())
        self.assertFalse(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

    def test_intron(self) -> None:
        v = VariantPosition("122-6")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (122, None, -6, None),
        )
        self.assertFalse(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

        v = VariantPosition("78+10")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr), (78, None, 10, None)
        )
        self.assertFalse(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

    def test_utr_intron(self) -> None:
        v = VariantPosition("*89+67")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr), (89, None, 67, True)
        )
        self.assertTrue(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

        v = VariantPosition("-127+6")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-127, None, 6, True),
        )
        self.assertTrue(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

        v = VariantPosition("*73-105")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (73, None, -105, True),
        )
        self.assertTrue(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())

        v = VariantPosition("-45-1")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-45, None, -1, True),
        )
        self.assertTrue(v.is_utr())
        self.assertTrue(v.is_intronic())
        self.assertFalse(v.is_protein())
        self.assertTrue(v.is_extended())


class TestObjectRepresentation(unittest.TestCase):
    def test_repr(self) -> None:
        position_strings = (
            "8",
            "92380",
            "*8",
            "-80",
            "122-6",
            "78+10",
            "*89+67",
            "-127+6",
            "*73-105",
            "-45-1",
            "Cys234",
            "Ala9",
        )
        for s in position_strings:
            with self.subTest(s=s):
                v = VariantPosition(s)
                self.assertEqual(s, repr(v))


# TODO: add amino acid variants
class TestComparisons(unittest.TestCase):
    def setUp(self) -> None:
        sorted_position_strings = (
            "-45-1",
            "-12",
            "8",
            "99",
            "99+88",
            "99+122",
            "100-12",
            "100",
            "101",
            "202-12",
            "202-1",
            "202",
            "*1",
            "*73-105",
        )

        self.sorted_variants = [VariantPosition(p) for p in sorted_position_strings]

        # pairwise itertools recipe
        a, b = itertools.tee(self.sorted_variants)
        next(b, None)
        self.sorted_variant_pairs = zip(a, b)

    def test_eq(self) -> None:
        for v in self.sorted_variants:
            with self.subTest(v=v):
                self.assertEqual(v, v)

    def test_ne(self) -> None:
        for v1, v2 in self.sorted_variant_pairs:
            with self.subTest(v1=v1, v2=v2):
                self.assertNotEqual(v1, v2)

    def test_lt(self) -> None:
        for v1, v2 in self.sorted_variant_pairs:
            with self.subTest(v1=v1, v2=v2):
                self.assertLess(v1, v2)

    def test_sorting(self) -> None:
        for _ in range(10):
            with self.subTest():
                shuffled_variants = self.sorted_variants.copy()
                while shuffled_variants == self.sorted_variants:
                    random.shuffle(shuffled_variants)
                self.assertListEqual(self.sorted_variants, sorted(shuffled_variants))


# TODO: add amino acid variants
class TestAdjacency(unittest.TestCase):
    def test_adjacent_pairs(self) -> None:
        adjacent_pairs = (
            ("-45-2", "-45-1"),
            ("-45-1", "-45"),
            ("-12", "-13"),
            ("-1", "1"),
            ("8", "9"),
            ("202-1", "202"),
            ("99", "99+1"),
            ("99+88", "99+89"),
            ("100-12", "100-11"),
            ("100", "101"),
            ("*1", "*2"),
            ("*73-1", "*73"),
        )
        for s1, s2 in adjacent_pairs:
            v1 = VariantPosition(s1)
            v2 = VariantPosition(s2)
            with self.subTest(v1=v1, v2=v2):
                self.assertTrue(v1.is_adjacent(v2))
            with self.subTest(v1=v1, v2=v2):
                self.assertTrue(v2.is_adjacent(v1))

    def test_not_adjacent_to_self(self) -> None:
        position_strings = (
            "-45-1",
            "-12",
            "8",
            "99",
            "99+88",
            "99+122",
            "100-12",
            "100",
            "103",
            "202-12",
            "202-1",
            "205",
            "*1",
            "*12",
            "*73-105",
        )
        variants = [VariantPosition(s) for s in position_strings]
        for v in variants:
            with self.subTest(v=v):
                self.assertFalse(v.is_adjacent(v))

    def test_non_adjacent_pairs(self) -> None:
        position_strings = (
            "-45-1",
            "-12",
            "8",
            "99",
            "99+88",
            "99+122",
            "100-12",
            "103",
            "202-12",
            "202-1",
            "205",
            "*1",
            "*12",
            "*73-105",
        )
        variants = [VariantPosition(s) for s in position_strings]

        for v1, v2 in itertools.permutations(variants, 2):
            with self.subTest(v1=v1, v2=v2):
                self.assertFalse(v1.is_adjacent(v2))


class TestHypothesisRoundTrip(unittest.TestCase):
    """Any string generated from the position grammar should parse successfully
    and format back to exactly the input string."""

    @given(s=valid_position_strings())
    def test_round_trip(self, s: str) -> None:
        v = VariantPosition(s)
        self.assertEqual(s, repr(v))

    @given(s=utr_intron_positions())
    def test_nucleotide_positions_are_never_protein(self, s: str) -> None:
        v = VariantPosition(s)
        self.assertFalse(v.is_protein())

    @given(s=amino_acid_position_strings())
    def test_amino_acid_positions_are_never_extended(self, s: str) -> None:
        v = VariantPosition(s)
        self.assertTrue(v.is_protein())
        self.assertFalse(v.is_extended())


class TestHypothesisFuzz(unittest.TestCase):
    """Any arbitrary string should either be rejected with MaveHgvsParseError or
    parse into an object whose repr reproduces the input; no other exception
    should ever escape the constructor."""

    @given(
        s=st.text(
            alphabet=string.ascii_letters + string.digits + "*+-",
            max_size=12,
        )
    )
    def test_fuzz_never_raises_unexpected_errors(self, s: str) -> None:
        try:
            v = VariantPosition(s)
        except MaveHgvsParseError:
            return
        self.assertEqual(s, repr(v))


class TestHypothesisComparisons(unittest.TestCase):
    @given(s=valid_position_strings())
    def test_equal_to_self(self, s: str) -> None:
        v = VariantPosition(s)
        self.assertEqual(v, v)
        self.assertFalse(v < v)

    @given(s1=valid_position_strings(), s2=valid_position_strings())
    def test_equality_is_symmetric(self, s1: str, s2: str) -> None:
        v1 = VariantPosition(s1)
        v2 = VariantPosition(s2)
        self.assertEqual(v1 == v2, v2 == v1)

    @given(s1=valid_position_strings(), s2=valid_position_strings())
    def test_ordering_is_consistent(self, s1: str, s2: str) -> None:
        v1 = VariantPosition(s1)
        v2 = VariantPosition(s2)
        # total_ordering trichotomy: exactly one of <, ==, > holds
        outcomes = (v1 < v2, v1 == v2, v1 > v2)
        self.assertEqual(1, sum(outcomes))
        # antisymmetry between the two directions
        self.assertEqual(v1 < v2, v2 > v1)


class TestHypothesisAdjacency(unittest.TestCase):
    @given(s=valid_position_strings())
    def test_never_adjacent_to_self(self, s: str) -> None:
        v = VariantPosition(s)
        self.assertFalse(v.is_adjacent(v))

    @given(s1=valid_position_strings(), s2=valid_position_strings())
    def test_adjacency_is_symmetric(self, s1: str, s2: str) -> None:
        v1 = VariantPosition(s1)
        v2 = VariantPosition(s2)
        self.assertEqual(v1.is_adjacent(v2), v2.is_adjacent(v1))

    @given(pair=adjacent_position_pairs())
    def test_generated_adjacent_pairs_are_adjacent(self, pair: tuple[str, str]) -> None:
        s1, s2 = pair
        v1 = VariantPosition(s1)
        v2 = VariantPosition(s2)
        self.assertTrue(v1.is_adjacent(v2))
        self.assertTrue(v2.is_adjacent(v1))


if __name__ == "__main__":
    unittest.main()
