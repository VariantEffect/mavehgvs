import unittest

from mavehgvs.util import parse_variant_strings
from mavehgvs.variant import Variant


class TestParseVariantStrings(unittest.TestCase):
    def test_sets_error_strings_for_invalid_items(self) -> None:
        invalid_variant_strings = [
            "g.Glu27Trp",
            "p.27Glu>Trp",
            "p.122-6T>A",
            "G>A",
            "22G>A",
            "G.44del",
            "a.78+5_78+10del",
            "77dup",
            "n.Pro12_Gly18dup",
            "g.22_23insauc",
            "g.25_24del",
            "g.25_24ins",
            "r.43-6_595+12delinsctt",
            "x.=",
            "c.(=)",
        ]

        for s in invalid_variant_strings:
            with self.subTest(s=s):
                valid, invalid = parse_variant_strings([s])
                self.assertIsNone(valid[0])
                self.assertIsInstance(invalid[0], str)

    def test_sets_variant_for_valid_items(self) -> None:
        valid_variant_strings = [
            "p.Glu27Trp",
            "c.122-6T>A",
            "g.44del",
            "c.78+5_78+10del",
            "c.77dup",
            "p.Pro12_Gly18dup",
            "p.Ala12_Pro13insGlyProCys",
            "r.22_23insauc",
            "c.43-6_595+12delinsCTT",
            "p.Ile71_Cys80delinsSer",
            "p.=",
            "c.=",
            "p.(=)",
        ]

        for s in valid_variant_strings:
            with self.subTest(s=s):
                valid, invalid = parse_variant_strings([s])
                self.assertIsInstance(valid[0], Variant)
                self.assertIsNone(invalid[0])

    def test_validates_with_targetseq(self) -> None:
        targetseq = "ACGT"
        valid_variant_strings = ["c.1A>T", "c.3G>C", "c.[1A>T;3G>C]"]
        invalid_variant_strings = ["c.1C>T", "c.3T>C", "c.[1A>T;3T>C]", "c.5A>G"]

        for s in valid_variant_strings:
            with self.subTest(s=s, targetseq=targetseq):
                valid, invalid = parse_variant_strings([s], targetseq=targetseq)
                self.assertIsInstance(valid[0], Variant)
                self.assertIsNone(invalid[0])

        for s in invalid_variant_strings:
            with self.subTest(s=s, targetseq=targetseq):
                valid, invalid = parse_variant_strings([s], targetseq=targetseq)
                self.assertIsNone(valid[0])
                self.assertIsInstance(invalid[0], str)

    def test_validates_expected_prefix(self) -> None:
        valid_variant_strings = ["p.Glu27Trp", "c.122-6T>A", "r.22_23insauc"]

        for s in valid_variant_strings:
            p = s[0]
            with self.subTest(s=s, p=p):
                valid, invalid = parse_variant_strings([s], expected_prefix=p)
                self.assertIsInstance(valid[0], Variant)
                self.assertIsNone(invalid[0])

        for s in valid_variant_strings:
            p = "g"
            with self.subTest(s=s, p=p):
                valid, invalid = parse_variant_strings([s], expected_prefix=p)
                self.assertIsNone(valid[0])
                self.assertIsInstance(invalid[0], str)

    def test_valid_expected_prefixes_only(self) -> None:
        valid_prefixes = list("cgmnopr")
        invalid_prefixes = list("CGMNOPRx.4ab?")
        variant = "p.Glu27Trp"

        for p in valid_prefixes:
            with self.subTest(p=p):
                parse_variant_strings([variant], expected_prefix=p)

        for p in invalid_prefixes:
            with self.subTest(p=p):
                with self.assertRaises(ValueError):
                    parse_variant_strings([variant], expected_prefix=p)
