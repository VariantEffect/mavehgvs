import unittest
from mavehgvs.variant import Variant


class TestCreateSingleVariantFromString(unittest.TestCase):
    def test_sub(self) -> None:
        variant_strings = [
            "p.Glu27Trp",
            "p.Ter345Lys",
            "p.Cys22=",
            "g.48C>A",
            "c.122-6T>A",
            "c.*33G>C",
            "r.22g>u",
            "r.33+12a>c",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_del(self) -> None:
        variant_strings = [
            "g.44del",
            "c.78+5_78+10del",
            "c.1_95del",
            "p.Gly18del",
            "p.Gln7_Asn19del",
            "r.34_36del",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_dup(self) -> None:
        variant_strings = [
            "g.22_24dup",
            "c.77dup",
            "c.101+1_101+7dup",
            "p.Pro12_Gly18dup",
            "p.Cys5dup",
            "r.12dup",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_ins(self) -> None:
        variant_strings = [
            "g.234_235insT",
            "c.84_85insCTG",
            "c.99+6_99+7insA",
            "p.His7_Gln8insSer",
            "p.Ala12_Pro13insGlyProCys",
            "r.22_23insauc",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_delins(self) -> None:
        variant_strings = [
            "g.22delinsAACG",
            "c.83_85delinsT",
            "c.43-6_595+12delinsCTT",
            "p.Ile71_Cys80delinsSer",
            "p.His44delinsValProGlyGlu",
            "r.92delinsgac",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_target_identical(self) -> None:
        variant_strings = [f"{prefix}.=" for prefix in "gmo" "cn" "r"]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_target_identical())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))


class TestCreateMultiVariantFromString(unittest.TestCase):
    pass


class TestCreateSingleVariantFromValues(unittest.TestCase):
    pass


class TestCreateMultiVariantFromValues(unittest.TestCase):
    pass


class TestTargetSequenceValidation(unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
