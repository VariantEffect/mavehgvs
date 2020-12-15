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
