import unittest

from mavehgvs.exceptions import MaveHgvsParseError
from mavehgvs.variant import Variant
from mavehgvs.position import VariantPosition


class TestCreateSingleVariantFromString(unittest.TestCase):
    def test_invalid_raises_error(self) -> None:
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
            "c.1_3=",
            "c.12=",
            "g.88_99=",
            "c.43-6_595+12=",
            "p.Glu27fs",
        ]

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
            "p.Pro12_Gly18insGlyProAla",
            "g.22_23insauc",
            "g.25_24del",
            "g.25_24ins",
            "r.22_24insauc",
            "r.43-6_595+12delinsctt",
            "x.=",
            "c.(=)",
            "p.(Gly24=)",
            "p.Gly24(=)",
            "p.Arg12LysfsTer18",
            "p.Glu27fs*?",
        ]

        for s in valid_variant_strings:
            with self.subTest(s=s):
                Variant(s)  # should pass

        for s in invalid_variant_strings:
            with self.subTest(s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s)

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
            "p.=",
            "p.(=)",
            "n.=",
            "c.1_3=",
            "c.12=",
            "g.88_99=",
            "c.43-6_595+12=",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_fs(self) -> None:
        variant_strings = ["p.Glu27fs"]

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
                self.assertEqual(s, str(v))

    def test_target_identical(self) -> None:
        identical_variant_strings = [
            *[f"{prefix}.=" for prefix in tuple("gmo" "cn" "r")],
            "p.(=)",
            "c.1_3=",
        ]

        non_identical_variant_strings = [
            "p.Ter345Lys",
            "p.Cys22=",
            "g.48C>A",
            "c.122-6T>A",
            "g.22delinsAACG",
            "c.83_85delinsT",
        ]

        for s in identical_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_target_identical())

        for s in non_identical_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertFalse(v.is_target_identical())

    def test_synonymous(self) -> None:
        synonymous_variant_strings = ["p.Gly24=", "p.=", "p.(=)"]

        nonsynonymous_variant_strings = ["p.Ter345Lys", "c.=", "g.48C>A"]

        for s in synonymous_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_synonymous())

        for s in nonsynonymous_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertFalse(v.is_synonymous())

    def test_relaxed_ordering(self):
        variant_tuples = [
            ("c.78+10_78+5del", "c.78+5_78+10del"),
            ("c.80_77dup", "c.77_80dup"),
            ("p.Gly18_Pro12dup", "p.Pro12_Gly18dup"),
            ("p.Pro13_Ala12insGlyProCys", "p.Ala12_Pro13insGlyProCys"),
            ("r.23_22insauc", "r.22_23insauc"),
            ("c.595+12_43-6delinsCTT", "c.43-6_595+12delinsCTT"),
            ("p.Cys80_Ile71delinsSer", "p.Ile71_Cys80delinsSer"),
            ("c.3_1=", "c.1_3="),
            ("g.99_88=", "g.88_99="),
            ("c.595+12_43-6=", "c.43-6_595+12="),
        ]

        for v, s in variant_tuples:
            with self.subTest(v=v, s=s):
                self.assertEqual(str(Variant(v, relaxed_ordering=True)), s)


class TestCreateMultiVariantFromString(unittest.TestCase):
    def test_creation(self):
        variant_strings = [
            "p.[Glu27Trp;Ter345Lys]",
            "p.[Glu27Trp;Lys212fs]",
            "p.[Gly18del;Glu27Trp;Ter345Lys]",
            "p.[Gln7_Asn19del;Glu27Trp;Ter345Lys]",
            "c.[1_35del;78+5_78+10del;122T>A]",
        ]

        invalid_variant_strings = [
            "p.[Glu27Trp;=;Ter345Lys]",
            "p.[(=);Gly18del;Glu27Trp;Ter345Lys]",
            "c.[12T>A;=;78+5_78+10del]",
            "c.[1_3=;12T>A;78+5_78+10del]",
            "p.[Glu27fs;Arg48Lys]",
            "p.[Glu27fs;Arg48fs]",
        ]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

        for s in invalid_variant_strings:
            with self.subTest(s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s)

    def test_ordering(self):
        variant_string_tuples = [
            ("p.[Gly345Lys;Glu27Trp]", "p.[Glu27Trp;Gly345Lys]"),
            ("p.[Glu27Trp;Gly18del;Ter345Lys]", "p.[Gly18del;Glu27Trp;Ter345Lys]"),
            ("c.[122T>A;1_35del;78+5_78+10del]", "c.[1_35del;78+5_78+10del;122T>A]"),
        ]

        for s, _ in variant_string_tuples:
            with self.subTest(s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, relaxed_ordering=False)

        for s, s_ordered in variant_string_tuples:
            with self.subTest(s=s):
                # Should pass creation
                Variant(s, relaxed_ordering=True)

        for s, s_ordered in variant_string_tuples:
            with self.subTest(s=s):
                v = Variant(s, relaxed_ordering=True)
                self.assertEqual(s_ordered, str(v))

    def test_overlaps(self):
        invalid_variant_strings = [
            "p.[Glu27Trp;Glu27Trp]",
            "p.[Glu27Trp;Glu27Tyr]",
            "p.[Pro27Trp;Glu27Tyr]",
            "p.[Gly18del;Gly18Tyr]",
            "p.[Gln7_Asn19del;Glu13Trp]",
            "p.[Glu13Trp;Gln7_Asn19del]",
            "p.[Gln7_Asn19del;Glu13Trp;Ter345Lys]",
            "c.[1_95del;78+5_78+10del;122T>A]",
            "c.[1_95del;22T>A]",
            "n.[22G>A;22G>T]",
        ]

        for s in invalid_variant_strings:
            with self.subTest(s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s)


class TestCreateSingleVariantFromValues(unittest.TestCase):
    def test_equal(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "equal",
                    "prefix": "p",
                    "position": "27",
                    "target": "Glu",
                },
                "p.Glu27=",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "c",
                    "start_position": "12",
                    "end_position": "12",
                },
                "c.12=",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "c",
                    "start_position": "1",
                    "end_position": "3",
                },
                "c.1_3=",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_sub(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "sub",
                    "prefix": "p",
                    "position": 27,
                    "target": "Glu",
                    "variant": "Trp",
                },
                "p.Glu27Trp",
            ),
            (
                {
                    "variant_type": "sub",
                    "prefix": "c",
                    "position": "122-6",
                    "target": "T",
                    "variant": "A",
                },
                "c.122-6T>A",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_fs(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "fs",
                    "prefix": "p",
                    "position": 27,
                    "target": "Glu",
                },
                "p.Glu27fs",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_ins(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "ins",
                    "prefix": "p",
                    "start_position": 12,
                    "start_target": "Ala",
                    "end_position": 13,
                    "end_target": "Pro",
                    "variant": "GlyProCys",
                },
                "p.Ala12_Pro13insGlyProCys",
            ),
            (
                {
                    "variant_type": "ins",
                    "prefix": "r",
                    "start_position": 22,
                    "end_position": 23,
                    "variant": "auc",
                },
                "r.22_23insauc",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_del(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "del",
                    "prefix": "g",
                    "start_position": 44,
                    "end_position": 44,
                },
                "g.44del",
            ),
            (
                {
                    "variant_type": "del",
                    "prefix": "c",
                    "start_position": "78+5",
                    "end_position": "78+10",
                },
                "c.78+5_78+10del",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_dup(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "dup",
                    "prefix": "c",
                    "start_position": 77,
                    "end_position": 77,
                },
                "c.77dup",
            ),
            (
                {
                    "variant_type": "dup",
                    "prefix": "p",
                    "start_position": 12,
                    "start_target": "Pro",
                    "end_position": 18,
                    "end_target": "Gly",
                },
                "p.Pro12_Gly18dup",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_delins(self):
        valid_dict_tuples = [
            (
                {
                    "variant_type": "delins",
                    "prefix": "c",
                    "start_position": "43-6",
                    "end_position": "595+12",
                    "variant": "CTT",
                },
                "c.43-6_595+12delinsCTT",
            ),
            (
                {
                    "variant_type": "delins",
                    "prefix": "p",
                    "start_position": 71,
                    "start_target": "Ile",
                    "end_position": 80,
                    "end_target": "Cys",
                    "variant": "Ser",
                },
                "p.Ile71_Cys80delinsSer",
            ),
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

    def test_extra_keys(self):
        invalid_dicts = [
            {
                "variant_type": "sub",
                "prefix": "p",
                "position": 27,
                "target": "Glu",
                "variant": "Trp",
                "bonus": "data",
            },
            {
                "variant_type": "sub",
                "prefix": "c",
                "position": "122-6",
                "start_target": "T",
                "target": "T",
                "variant": "A",
            },
            {
                "variant_type": "delins",
                "prefix": "p",
                "start_target": "Ile",
                "end_position": 80,
                "end_target": "Cys",
                "variant": "Ser",
                "position": "Ala",
            },
        ]

        for d in invalid_dicts:
            with self.subTest(d=d):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(d)

    def test_missing_keys(self):
        invalid_dicts = [
            {"prefix": "p", "position": 27, "target": "Glu", "variant": "Trp"},
            {"variant_type": "sub", "position": "122-6", "target": "T", "variant": "A"},
            {
                "variant_type": "delins",
                "prefix": "p",
                "start_target": "Ile",
                "end_position": 80,
                "end_target": "Cys",
                "variant": "Ser",
            },
        ]

        for d in invalid_dicts:
            with self.subTest(d=d):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(d)

    def test_invalid_keys(self):
        invalid_dicts = [
            {
                "variant_type": "equal",
                "prefix": "p",
                "start_position": "27",
                "end_position": "27",
                "target": "Glu",
            },
            {"variant_type": "dup", "prefix": "c", "position": 77},
            {
                "variant_type": "test",
                "prefix": "c",
                "start_position": 77,
                "end_position": 77,
            },
            {
                "variant_type": "fs",
                "prefix": "c",
                "position": "12",
                "target": "T",
            },
        ]

        for d in invalid_dicts:
            with self.subTest(d=d):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(d)

    def test_invalid_type(self):
        invalid_values = [1234, None, 5.55, ("p", "Ile", 80, "Cys")]

        for v in invalid_values:
            with self.subTest(v=v):
                with self.assertRaises(ValueError):
                    Variant(v)


class TestCreateMultiVariantFromValues(unittest.TestCase):
    def test_create_multivariant(self):
        valid_dict_tuples = [
            (
                [
                    {
                        "variant_type": "sub",
                        "prefix": "p",
                        "position": 27,
                        "target": "Glu",
                        "variant": "Trp",
                    },
                    {
                        "variant_type": "delins",
                        "prefix": "p",
                        "start_position": 71,
                        "start_target": "Ile",
                        "end_position": 80,
                        "end_target": "Cys",
                        "variant": "Ser",
                    },
                ],
                "p.[Glu27Trp;Ile71_Cys80delinsSer]",
            ),
            (
                [
                    {
                        "variant_type": "dup",
                        "prefix": "c",
                        "start_position": 77,
                        "end_position": 77,
                    },
                    {
                        "variant_type": "sub",
                        "prefix": "c",
                        "position": "122-6",
                        "target": "T",
                        "variant": "A",
                    },
                ],
                "c.[77dup;122-6T>A]",
            ),
        ]

        invalid_dicts = [
            [
                {
                    "variant_type": "sub",
                    "position": 27,
                    "target": "Glu",
                    "variant": "Trp",
                },
                {
                    "variant_type": "delins",
                    "prefix": "p",
                    "start_position": 71,
                    "start_target": "Ile",
                    "end_position": 80,
                    "end_target": "Cys",
                    "variant": "Ser",
                },
            ],
            [
                {
                    "variant_type": "sub",
                    "prefix": "p",
                    "position": 27,
                    "target": "Glu",
                    "variant": "Trp",
                },
                {
                    "variant_type": "sub",
                    "prefix": "c",
                    "position": "122-6",
                    "target": "T",
                    "variant": "A",
                },
            ],
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

        for d in invalid_dicts:
            with self.subTest(d=d):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(d)


class TestTargetSequenceValidation(unittest.TestCase):
    def test_valid_dna_equal(self):
        variant_tuples = [("ACGT", "c.1_2="), ("ACGT", "c.4="), ("ACGT", "c.=")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_dna_equal(self):
        variant_tuples = [("ACGT", "c.4_5="), ("ACGT", "c.10=")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_matching_dna_substitution(self):
        variant_tuples = [
            ("ACGT", "c.1A>T"),
            ("ACGT", "c.3G>C"),
            ("ACGT", "c.[1A>T;3G>C]"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_nonmatching_dna_substitution(self):
        variant_tuples = [
            ("ACGT", "c.1C>T"),
            ("ACGT", "c.3T>C"),
            ("ACGT", "c.[1A>T;3T>C]"),
            ("ACGT", "c.5A>G"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_valid_dna_del(self):
        variant_tuples = [("ACGT", "c.1_3del"), ("ACGT", "c.4del")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_dna_del(self):
        variant_tuples = [
            ("ACGT", "c.1_5del"),
            ("ACGT", "c.6_8del"),
            ("ACGT", "c.7del"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_valid_dna_dup(self):
        variant_tuples = [("ACGT", "c.1_3dup"), ("ACGT", "c.4dup")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_dna_dup(self):
        variant_tuples = [
            ("ACGT", "c.1_5dup"),
            ("ACGT", "c.6_8dup"),
            ("ACGT", "c.7dup"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_valid_dna_ins(self):
        variant_tuples = [("ACGT", "c.1_2insAAA"), ("ACGT", "c.3_4insT")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_dna_ins(self):
        variant_tuples = [("ACGT", "c.4_5insA"), ("ACGT", "c.10_11insTCG")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_valid_dna_delins(self):
        variant_tuples = [("ACGT", "c.1_2delinsA"), ("ACGT", "c.4delinsTAAGC")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_dna_delins(self):
        variant_tuples = [("ACGT", "c.4_5delinsA"), ("ACGT", "c.10_delinsTCG")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_valid_protein_equal(self):
        variant_tuples = [("RCQY", "p.Arg1="), ("RCQY", "p.Tyr4="), ("RCQY", "p.=")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_invalid_protein_equal(self):
        variant_tuples = [("RCQY", "p.Trp5=")]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_matching_protein_substitution(self):
        variant_tuples = [
            ("RCQY", "p.Arg1Ala"),
            ("RCQY", "p.Gln3Trp"),
            ("RCQY", "p.[Arg1Ala;Gln3Trp]"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_nonmatching_protein_substitution(self):
        variant_tuples = [
            ("RCQY", "p.Cys1Ala"),
            ("RCQY", "p.Ala3Trp"),
            ("RCQY", "p.[Arg1Ala;Cys3Trp]"),
            ("RCQY", "p.Asp5Glu"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_matching_protein_fs(self):
        variant_tuples = [
            ("RCQY", "p.Arg1fs"),
            ("RCQY", "p.Gln3fs"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_nonmatching_protein_fs(self):
        variant_tuples = [
            ("RCQY", "p.Cys1fs"),
            ("RCQY", "p.Ala3fs"),
            ("RCQY", "p.Asp5fs"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_matching_protein_indel(self):
        variant_tuples = [
            ("RCQY", "p.Arg1del"),
            ("RCQY", "p.Arg1_Gln3dup"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))

    def test_nonmatching_protein_indel(self):
        variant_tuples = [
            ("RCQY", "p.Cys1del"),
            ("RCQY", "p.Arg1_Asp3dup"),
            ("RCQY", "p.Asp5del"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(s, targetseq=target)

    def test_skips_extended(self):
        variant_tuples = [
            ("ACGT", "c.1+3A>T"),
            ("ACGT", "c.*33G>C"),
            ("ACGT", "c.43-6_595+12delinsCTT"),
        ]

        for target, s in variant_tuples:
            with self.subTest(target=target, s=s):
                v = Variant(s, targetseq=target)
                self.assertEqual(s, str(v))


class TestMiscMethods(unittest.TestCase):
    def test_is_multi_variant(self):
        single_variant_strings = [
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
        ]

        multi_variant_strings = []

        for s in single_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertFalse(v.is_multi_variant())

        for s in multi_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_multi_variant())

    def test_uses_extended_positions(self):
        non_extended_variant_strings = [
            "p.Glu27Trp",
            "g.44del",
            "c.77dup",
            "p.Pro12_Gly18dup",
            "p.Ala12_Pro13insGlyProCys",
            "r.22_23insauc",
            "r.22g>u",
            "p.Ile71_Cys80delinsSer",
            "p.=",
            "p.[Pro12_Gly18dup;Glu27Trp]",
            "r.[22g>u;35del]",
        ]

        extended_variant_strings = [
            "c.122-6T>A",
            "c.78+5_78+10del",
            "c.43-6_595+12delinsCTT",
            "c.*33G>C",
            "r.33+12a>c",
            "c.[12G>T;122-6T>A]",
            "c.[43-6_595+12delinsCTT;*33G>C]",
        ]

        for s in non_extended_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertFalse(v.uses_extended_positions())

        for s in extended_variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.uses_extended_positions())


# TODO: multi-variant test cases
class TestMiscProperties(unittest.TestCase):
    def test_prefix(self):
        variant_tuples = [(prefix, f"{prefix}.=") for prefix in tuple("gmo" "cn" "r")]

        for p, s in variant_tuples:
            with self.subTest(p=p, s=s):
                v = Variant(s)
                self.assertEqual(p, v.prefix)

    def test_variant_type(self):
        variant_tuples = [
            ("sub", "p.Glu27Trp"),
            ("sub", "c.122-6T>A"),
            ("fs", "p.Glu27fs"),
            ("del", "g.44del"),
            ("del", "c.78+5_78+10del"),
            ("dup", "c.77dup"),
            ("dup", "p.Pro12_Gly18dup"),
            ("ins", "p.Ala12_Pro13insGlyProCys"),
            ("ins", "r.22_23insauc"),
            ("delins", "c.43-6_595+12delinsCTT"),
            ("delins", "p.Ile71_Cys80delinsSer"),
        ]

        for t, s in variant_tuples:
            with self.subTest(t=t, s=s):
                v = Variant(s)
                self.assertEqual(t, v.variant_type)

    def test_position(self):
        variant_tuples = [
            (VariantPosition("Glu27"), "p.Glu27Trp"),
            (VariantPosition("Glu27"), "p.Glu27fs"),
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

    def test_sequence(self):
        variant_tuples = [
            (("Glu", "Trp"), "p.Glu27Trp"),
            (("T", "A"), "c.122-6T>A"),
            (None, "p.Glu27fs"),
            (None, "g.44del"),
            (None, "c.78+5_78+10del"),
            (None, "c.77dup"),
            (None, "p.Pro12_Gly18dup"),
            ("GlyProCys", "p.Ala12_Pro13insGlyProCys"),
            ("auc", "r.22_23insauc"),
            ("CTT", "c.43-6_595+12delinsCTT"),
            ("Ser", "p.Ile71_Cys80delinsSer"),
        ]

        for seq, s in variant_tuples:
            with self.subTest(seq=seq, s=s):
                v = Variant(s)
                self.assertEqual(seq, v.sequence)

    def test_target_id(self):
        variant_tuples = [
            (None, "p.Glu27Trp"),
            (None, "c.122-6T>A"),
            ("GeneX", "GeneX:p.Glu27Trp"),
            ("YFG1", "YFG1:c.122-6T>A"),
            ("ENST00000471181.7", "ENST00000471181.7:c.122-6T>A"),
            ("NM_007294.4", "NM_007294.4:c.122-6T>A"),
        ]

        for t, s in variant_tuples:
            with self.subTest(t=t, s=s):
                v = Variant(s)
                self.assertEqual(t, v.target_id)

        for _, s in variant_tuples:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))


if __name__ == "__main__":
    unittest.main()
