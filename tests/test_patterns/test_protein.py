import unittest
import re
from mavehgvs.patterns.protein import (
    pro_equal,
    pro_sub,
    pro_fs,
    pro_del,
    pro_dup,
    pro_ins,
    pro_delins,
    pro_variant,
    pro_single_variant,
    pro_multi_variant,
)
from . import build_multi_variants


class TestProteinEqual(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_equal, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "(=)",
            "Cys22=",
        ]

        cls.invalid_strings = ["=22", "Arg18(=)", "Cys-22", "=="]

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


class TestProteinSub(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_sub, flags=re.ASCII)

        cls.valid_strings = ["Glu27Trp", "Ter345Lys"]

        cls.invalid_strings = [
            "22A>T",
            "Xaa12Arg",
            "Arg21Xaa",
            "Pro17*",
            "*345Lys",
            "(Glu27Trp)",
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


class TestProteinFs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_fs, flags=re.ASCII)

        cls.valid_strings = ["Glu27fs"]

        cls.invalid_strings = [
            "=fs",
            "Arg12LysfsTer18",
            "Arg12Lysfs*18",
            "Glu27fs*?",
            "(Glu27fs)",
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


class TestProteinDel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_del, flags=re.ASCII)

        cls.valid_strings = [
            "Gly18del",
            "Gln7_Asn19del",
        ]

        cls.invalid_strings = ["=del", "18del", "122_128del", "(Gly18del)"]

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


class TestProteinDup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_dup, flags=re.ASCII)

        cls.valid_strings = [
            "Cys5dup",
            "Pro12_Gly18dup",
        ]

        cls.invalid_strings = ["=dup", "18dup", "122_128dup", "(Cys5dup)"]

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


class TestProteinIns(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_ins, flags=re.ASCII)

        cls.valid_strings = [
            "His7_Gln8insSer",
            "Ala12_Pro13insGlyProCys",
        ]

        cls.invalid_strings = [
            "(His7_Gln8insSer)",
            "(His7_Gln8insX)",
            "(Ala12_Pro13ins(2))",
            "His7_Gln8ins?",
            "His7_Gln8insXaa",
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


class TestProteinDelins(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_delins, flags=re.ASCII)

        cls.valid_strings = [
            "Ile71_Cys80delinsSer",
            "His44delinsValProGlyGlu",
        ]

        cls.invalid_strings = ["(Ile71_Cys80delinsSer)", "Ile71_Cys80delinsXaa"]

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


class TestProteinVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_variant, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "(=)",
            "Cys22=",
            "Glu27Trp",
            "Ter345Lys",
            "Glu27fs",
            "Gly18del",
            "Gln7_Asn19del",
            "Cys5dup",
            "Pro12_Gly18dup",
            "His7_Gln8insSer",
            "Ala12_Pro13insGlyProCys",
            "Ile71_Cys80delinsSer",
            "His44delinsValProGlyGlu",
        ]

        cls.invalid_strings = [
            "=22",
            "Arg18(=)",
            "Cys-22",
            "==",
            "22A>T",
            "Xaa12Arg",
            "Arg21Xaa",
            "Pro17*",
            "*345Lys",
            "(Glu27Trp)",
            "=fs",
            "Arg12LysfsTer18",
            "Arg12Lysfs*18",
            "Glu27fs*?",
            "(Glu27fs)",
            "=del",
            "18del",
            "122_128del",
            "(Gly18del)",
            "=dup",
            "18dup",
            "122_128dup",
            "(Cys5dup)",
            "(His7_Gln8insSer)",
            "(His7_Gln8insX)",
            "(Ala12_Pro13ins(2))",
            "His7_Gln8ins?",
            "His7_Gln8insXaa",
            "(Ile71_Cys80delinsSer)",
            "Ile71_Cys80delinsXaa",
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


class TestProteinSingleVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_single_variant, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "(=)",
            "Cys22=",
            "Glu27Trp",
            "Ter345Lys",
            "Glu27fs",
            "Gly18del",
            "Gln7_Asn19del",
            "Cys5dup",
            "Pro12_Gly18dup",
            "His7_Gln8insSer",
            "Ala12_Pro13insGlyProCys",
            "Ile71_Cys80delinsSer",
            "His44delinsValProGlyGlu",
        ]

        cls.invalid_strings = [
            "=22",
            "Arg18(=)",
            "Cys-22",
            "==",
            "22A>T",
            "Xaa12Arg",
            "Arg21Xaa",
            "Pro17*",
            "*345Lys",
            "(Glu27Trp)",
            "=fs",
            "Arg12LysfsTer18",
            "Arg12Lysfs*18",
            "Glu27fs*?",
            "(Glu27fs)",
            "=del",
            "18del",
            "122_128del",
            "(Gly18del)",
            "=dup",
            "18dup",
            "122_128dup",
            "(Cys5dup)",
            "(His7_Gln8insSer)",
            "(His7_Gln8insX)",
            "(Ala12_Pro13ins(2))",
            "His7_Gln8ins?",
            "His7_Gln8insXaa",
            "(Ile71_Cys80delinsSer)",
            "Ile71_Cys80delinsXaa",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                v = f"p.{s}"
                self.assertIsNotNone(
                    self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                v = f"p.{s}"
                self.assertIsNone(
                    self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                )


class TestProteinMultiVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(pro_multi_variant, flags=re.ASCII)

        single_valid_strings = [
            "=",
            "(=)",
            "Cys22=",
            "Glu27Trp",
            "Ter345Lys",
            "Glu27fs",
            "Gly18del",
            "Gln7_Asn19del",
            "Cys5dup",
            "Pro12_Gly18dup",
            "His7_Gln8insSer",
            "Ala12_Pro13insGlyProCys",
            "Ile71_Cys80delinsSer",
            "His44delinsValProGlyGlu",
        ]

        single_invalid_strings = [
            "=22",
            "Arg18(=)",
            "Cys-22",
            "==",
            "22A>T",
            "Xaa12Arg",
            "Arg21Xaa",
            "Pro17*",
            "*345Lys",
            "(Glu27Trp)",
            "=fs",
            "Arg12LysfsTer18",
            "Arg12Lysfs*18",
            "Glu27fs*?",
            "(Glu27fs)",
            "=del",
            "18del",
            "122_128del",
            "(Gly18del)",
            "=dup",
            "18dup",
            "122_128dup",
            "(Cys5dup)",
            "(His7_Gln8insSer)",
            "(His7_Gln8insX)",
            "(Ala12_Pro13ins(2))",
            "His7_Gln8ins?",
            "His7_Gln8insXaa",
            "(Ile71_Cys80delinsSer)",
            "Ile71_Cys80delinsXaa",
        ]

        cls.valid_strings, cls.invalid_strings = build_multi_variants(
            single_valid_strings, single_invalid_strings
        )

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                v = f"p.[{s}]"
                self.assertIsNotNone(
                    self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                v = f"p.[{s}]"
                self.assertIsNone(
                    self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                )


if __name__ == "__main__":
    unittest.main()
