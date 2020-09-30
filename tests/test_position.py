import unittest
import itertools
import random
from mavehgvs.variant import VariantPosition


class TestObjectCreation(unittest.TestCase):
    def test_position_only(self) -> None:
        v = VariantPosition("8")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (8, None, None))

        v = VariantPosition("92380")
        self.assertTupleEqual(
            (v.position, v.intronic_position, v.utr), (92380, None, None)
        )

    def test_invalid_strings(self) -> None:
        variant_strings = (
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
        )
        for s in variant_strings:
            with self.subTest(s=s):
                with self.assertRaises(ValueError):
                    VariantPosition(s)

    def test_utr(self) -> None:
        v = VariantPosition("*8")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (8, None, True))

        v = VariantPosition("-80")
        self.assertTupleEqual(
            (v.position, v.intronic_position, v.utr), (-80, None, True)
        )

    def test_intron(self) -> None:
        v = VariantPosition("122-6")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (122, -6, None))

        v = VariantPosition("78+10")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (78, 10, None))

    def test_utr_intron(self) -> None:
        v = VariantPosition("*89+67")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (89, 67, True))

        v = VariantPosition("-127+6")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (-127, 6, True))

        v = VariantPosition("*73-105")
        self.assertTupleEqual(
            (v.position, v.intronic_position, v.utr), (73, -105, True)
        )

        v = VariantPosition("-45-1")
        self.assertTupleEqual((v.position, v.intronic_position, v.utr), (-45, -1, True))


class TestObjectRepresentation(unittest.TestCase):
    def test_repr(self) -> None:
        variant_strings = (
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
        )
        for s in variant_strings:
            with self.subTest(s=s):
                v = VariantPosition(s)
                self.assertEqual(s, repr(v))


class TestComparisons(unittest.TestCase):
    def setUp(self) -> None:
        sorted_variant_strings = (
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

        self.sorted_variants = [VariantPosition(p) for p in sorted_variant_strings]

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
        variant_strings = (
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
        variants = [VariantPosition(s) for s in variant_strings]
        for v in variants:
            with self.subTest(v=v):
                self.assertFalse(v.is_adjacent(v))

    def test_non_adjacent_pairs(self) -> None:
        variant_strings = (
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
        variants = [VariantPosition(s) for s in variant_strings]

        for v1, v2 in itertools.permutations(variants, 2):
            with self.subTest(v1=v1, v2=v2):
                self.assertFalse(v1.is_adjacent(v2))


if __name__ == "__main__":
    unittest.main()
