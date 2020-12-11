import unittest
import itertools
import random
from mavehgvs.position import VariantPosition
from mavehgvs.exceptions import MaveHgvsParseError


class TestObjectCreation(unittest.TestCase):
    def test_position_only(self) -> None:
        v = VariantPosition("8")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (8, None, None, None),
        )

        v = VariantPosition("92380")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (92380, None, None, None),
        )

    def test_amino_acid(self) -> None:
        v = VariantPosition("Gly8")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (8, "Gly", None, None),
        )

        v = VariantPosition("Cys92380")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (92380, "Cys", None, None),
        )

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

        v = VariantPosition("-80")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-80, None, None, True),
        )

    def test_intron(self) -> None:
        v = VariantPosition("122-6")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (122, None, -6, None),
        )

        v = VariantPosition("78+10")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr), (78, None, 10, None)
        )

    def test_utr_intron(self) -> None:
        v = VariantPosition("*89+67")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr), (89, None, 67, True)
        )

        v = VariantPosition("-127+6")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-127, None, 6, True),
        )

        v = VariantPosition("*73-105")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (73, None, -105, True),
        )

        v = VariantPosition("-45-1")
        self.assertTupleEqual(
            (v.position, v.amino_acid, v.intronic_position, v.utr),
            (-45, None, -1, True),
        )


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


if __name__ == "__main__":
    unittest.main()
