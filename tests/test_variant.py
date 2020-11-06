import unittest
from mavehgvs.variant import Variant
from mavehgvs.position import VariantPosition


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


class TestMiscMethods(unittest.TestCase):
    @unittest.expectedFailure
    def test_is_valid(self):
        self.assertTrue(False)

    @unittest.expectedFailure
    def test_is_multi_variant(self):
        self.assertTrue(False)

    @unittest.expectedFailure
    def test_uses_extended_positions(self):
        self.assertTrue(False)


# TODO: multi-variant test cases
class TestMiscProperties(unittest.TestCase):
    def test_prefix(self):
        variant_tuples = [(prefix, f"{prefix}.=") for prefix in "gmo" "cn" "r"]
        variant_tuples.append((None, "x.="))  # invalid variant

        for p, s in variant_tuples:
            with self.subTest(p=p, s=s):
                v = Variant(s)
                self.assertEqual(p, v.prefix)

    def test_variant_type(self):
        variant_tuples = [
            ("sub", "p.Glu27Trp"),
            ("sub", "c.122-6T>A"),
            ("del", "g.44del"),
            ("del", "c.78+5_78+10del"),
            ("dup", "c.77dup"),
            ("dup", "p.Pro12_Gly18dup"),
            ("ins", "p.Ala12_Pro13insGlyProCys"),
            ("ins", "r.22_23insauc"),
            ("delins", "c.43-6_595+12delinsCTT"),
            ("delins", "p.Ile71_Cys80delinsSer"),
            (None, "x.="),
        ]

        for t, s in variant_tuples:
            with self.subTest(t=t, s=s):
                v = Variant(s)
                self.assertEqual(t, v.variant_type)

    def test_position(self):
        variant_tuples = [
            (VariantPosition("Glu27"), "p.Glu27Trp"),
            (VariantPosition("122-6"), "c.122-6T>A"),
            (VariantPosition("44"), "g.44del"),
            ((VariantPosition("78+5"), VariantPosition("78+10")), "c.78+5_78+10del"),
            (VariantPosition("77"), "c.77dup"),
            ((VariantPosition("Pro12"), VariantPosition("Gly18")), "p.Pro12_Gly18dup"),
            (
                (VariantPosition("Ala12"), VariantPosition("Pro13")),
                "p.Ala12_Pro13insGlyProCys",
            ),
            ((VariantPosition("22"), VariantPosition("23")), "r.22_23insauc"),
            (
                (VariantPosition("43-6"), VariantPosition("595+12")),
                "c.43-6_595+12delinsCTT",
            ),
            (
                (VariantPosition("Ile71"), VariantPosition("Cys80")),
                "p.Ile71_Cys80delinsSer",
            ),
            (None, "x.="),
        ]

        for p, s in variant_tuples:
            with self.subTest(p=p, s=s):
                v = Variant(s)
                if isinstance(p, list):  # multi-variant
                    self.assertEqual(len(p), len(v.positions))
                    for q, vp in zip(p, v.positions):
                        if isinstance(q, tuple):
                            self.assertTupleEqual(q, vp)
                        else:
                            self.assertEqual(q, vp)
                if isinstance(p, tuple):
                    self.assertTupleEqual(p, v.positions)
                else:
                    self.assertEqual(p, v.positions)

    @unittest.expectedFailure
    def test_sequence(self):
        self.assertTrue(False)

    @unittest.expectedFailure
    def test_target_id(self):
        self.assertTrue(False)


if __name__ == "__main__":
    unittest.main()
