import unittest
import re
from mavehgvs.patterns.rna import (
    rna_equal,
    rna_sub,
    rna_del,
    rna_dup,
    rna_ins,
    rna_delins,
    rna_variant,
    rna_single_variant,
    rna_multi_variant,
)
from . import build_multi_variants


class TestRnaEqual(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_equal, flags=re.ASCII)

        cls.valid_strings = [
            "=",
        ]

        cls.invalid_strings = ["=22", "(=)", "=="]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaSub(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_sub, flags=re.ASCII)

        cls.valid_strings = ["22g>u", "33+12a>c"]

        cls.invalid_strings = [
            "spl",
            "33+12A>G",
            "22g>t",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaDel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_del, flags=re.ASCII)

        cls.valid_strings = ["34_36del", "17del", "27_27+12del", "101+1_101+7del"]

        cls.invalid_strings = ["=del", "=/9_12del", "(155_185)del", "34_36"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaDup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_dup, flags=re.ASCII)

        cls.valid_strings = ["12dup", "2_24dup", "101+1_101+7dup", "12-24_12-12dup"]

        cls.invalid_strings = ["=dup", "(78+1_79-1)_(124+1_125-1)dup"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaIns(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_ins, flags=re.ASCII)

        cls.valid_strings = [
            "22_23insauc",
            "17_18insa",
        ]

        cls.invalid_strings = [
            "(27_30)insu",
            "74_74insnnn",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaDelins(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_delins, flags=re.ASCII)

        cls.valid_strings = ["92delinsgac", "12_17delinsc"]

        cls.invalid_strings = ["234_235ins(10)", "(122_125)insg"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_variant, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "22g>u",
            "33+12a>c",
            "34_36del",
            "17del",
            "12dup",
            "2_24dup",
            "101+1_101+7dup",
            "22_23insauc",
            "17_18insa",
            "92delinsgac",
            "12_17delinsc",
        ]

        cls.invalid_strings = [
            "=22",
            "(=)",
            "==",
            "spl",
            "33+12A>G",
            "22g>t",
            "=del",
            "=/9_12del",
            "(155_185)del",
            "=dup",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(27_30)insu",
            "74_74insnnn",
            "234_235ins(10)",
            "(122_125)insg",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestRnaSingleVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_single_variant, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "22g>u",
            "33+12a>c",
            "34_36del",
            "17del",
            "12dup",
            "2_24dup",
            "101+1_101+7dup",
            "22_23insauc",
            "17_18insa",
            "92delinsgac",
            "12_17delinsc",
        ]

        cls.invalid_strings = [
            "=22",
            "(=)",
            "==",
            "spl",
            "33+12A>G",
            "22g>t",
            "=del",
            "=/9_12del",
            "(155_185)del",
            "=dup",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(27_30)insu",
            "74_74insnnn",
            "234_235ins(10)",
            "(122_125)insg",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                v = f"r.{s}"
                self.assertIsNotNone(
                    self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                v = f"r.{s}"
                self.assertIsNone(
                    self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                )


class TestRnaMultiVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(rna_multi_variant, flags=re.ASCII)

        single_valid_strings = [
            "=",
            "22g>u",
            "33+12a>c",
            "34_36del",
            "17del",
            "12dup",
            "2_24dup",
            "101+1_101+7dup",
            "22_23insauc",
            "17_18insa",
            "92delinsgac",
            "12_17delinsc",
        ]

        single_invalid_strings = [
            "=22",
            "(=)",
            "==",
            "spl",
            "33+12A>G",
            "22g>t",
            "=del",
            "=/9_12del",
            "(155_185)del",
            "=dup",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(27_30)insu",
            "74_74insnnn",
            "234_235ins(10)",
            "(122_125)insg",
        ]
        cls.valid_strings, cls.invalid_strings = build_multi_variants(
            single_valid_strings, single_invalid_strings
        )

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                v = f"r.[{s}]"
                self.assertIsNotNone(
                    self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                v = f"r.[{s}]"
                self.assertIsNone(
                    self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                )


if __name__ == "__main__":
    unittest.main()
