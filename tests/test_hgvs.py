from unittest import TestCase

from hgvsp import multi_variant_re, single_variant_re, any_variant_re
from hgvsp import Level, Event, infer_level, infer_type, is_multi


class TestInferLevel(TestCase):
    def test_infers_protein(self):
        self.assertEqual(Level.PROTEIN, infer_level("p"))

    def test_infers_dna(self):
        self.assertEqual(Level.DNA, infer_level("c"))
        self.assertEqual(Level.DNA, infer_level("n"))
        self.assertEqual(Level.DNA, infer_level("g"))
        self.assertEqual(Level.DNA, infer_level("m"))

    def test_infers_rna(self):
        self.assertEqual(Level.RNA, infer_level("r"))

    def test_infers_none(self):
        self.assertEqual(None, infer_level("t"))

    def test_none_null_value(self):
        for v in [False, None, ""]:
            self.assertIsNone(infer_level(v))


class TestInferType(TestCase):
    def test_infers_substitution(self):
        self.assertEqual(Event.SUBSTITUTION, infer_type("c.1A>G"))
        self.assertEqual(Event.SUBSTITUTION, infer_type("r.1a>u"))
        self.assertEqual(Event.SUBSTITUTION, infer_type("p.Gly4Leu"))

    def test_infers_insertion(self):
        self.assertEqual(Event.INSERTION, infer_type("c.240_241insAGG"))
        self.assertEqual(Event.INSERTION, infer_type("r.2949_2950insaaa"))
        self.assertEqual(Event.INSERTION, infer_type("p.Arg78_Gly79ins23"))

    def test_infers_deletion(self):
        self.assertEqual(Event.DELETION, infer_type("c.(?_-1)_(*1_?)del"))
        self.assertEqual(Event.DELETION, infer_type("r.6_8del"))
        self.assertEqual(Event.DELETION, infer_type("p.Val7=/del"))

    def test_infers_delins(self):
        self.assertEqual(Event.DELINS, infer_type("c.6775delinsGA"))
        self.assertEqual(Event.DELINS, infer_type("r.?_?delinsc"))
        self.assertEqual(Event.DELINS, infer_type("p.C28_L29delinsTG"))

    def test_infers_frame_shift(self):
        self.assertEqual(Event.FRAME_SHIFT, infer_type("Glu5ValfsTer5"))
        self.assertEqual(Event.FRAME_SHIFT, infer_type("p.Ile327Argfs*?"))

    def test_none_null_value(self):
        for v in [False, None, ""]:
            self.assertIsNone(infer_type(v))


class TestIsMulti(TestCase):
    def test_detects_multi_variant(self):
        self.assertTrue(is_multi("c.[123A>G;19del]"))
        self.assertTrue(is_multi("n.[123A>G;19del]"))
        self.assertTrue(is_multi("g.[123A>G;19del]"))
        self.assertTrue(is_multi("m.[123A>G;19del]"))

        self.assertTrue(is_multi("r.[123a>g;19del]"))
        self.assertTrue(is_multi("r.[123a>g,19del]"))

        self.assertTrue(is_multi("p.[His4_Gln5insAla;Cys28fs;Cys28delinsVal]"))

    def test_single_variant_returns_false(self):
        self.assertFalse(is_multi("c.123A>G"))
        self.assertFalse(is_multi("r.19del"))
        self.assertFalse(is_multi("r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]"))
        self.assertFalse(is_multi("p.His4_Gln5insAla"))
        self.assertFalse(is_multi("p.(His4_Gln5insAla)"))

    def test_incomplete_multi_returns_false(self):
        self.assertFalse(is_multi("r.[123a>g;19del;]"))
        self.assertFalse(is_multi("r.[123a>g,19del,]"))
        self.assertFalse(is_multi("c.[123a>g]"))
        self.assertFalse(is_multi("c.[]"))
        self.assertFalse(is_multi("g.[]"))
        self.assertFalse(is_multi("n.[]"))
        self.assertFalse(is_multi("m.[]"))
        self.assertFalse(is_multi("r.[]"))
        self.assertFalse(is_multi("p.[]"))

    def test_mixed_multi_variant_returns_false(self):
        self.assertFalse(is_multi("c.[1A>G;Lys4Gly]"))


class TestMatchMulti(TestCase):
    def test_no_match_invalid_prefix(self):
        self.assertIsNone(multi_variant_re.fullmatch("f.1A>G"))

    def test_can_validate_protein_multi(self):
        self.assertIsNotNone(
            multi_variant_re.fullmatch("p.[His4_Gln5insAla;Cys28fs;Cys28delinsVal]")
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch("p.[(His4_Gln5insAla);(Cys28fs);Cys28delinsVal]")
        )

    def test_can_validate_dna_multi(self):
        self.assertIsNotNone(multi_variant_re.fullmatch("c.[123A>G;19del]"))
        self.assertIsNotNone(multi_variant_re.fullmatch("g.[123A>G;19del]"))
        self.assertIsNotNone(multi_variant_re.fullmatch("m.[123A>G;19del]"))
        self.assertIsNotNone(multi_variant_re.fullmatch("n.[123A>G;19del]"))

    def test_can_validate_rna_multi(self):
        self.assertIsNotNone(multi_variant_re.fullmatch("r.[123a>g;19del]"))
        self.assertIsNotNone(multi_variant_re.fullmatch("r.[123a>g,19del]"))

    def test_does_not_match_invalid_multi(self):
        self.assertIsNone(multi_variant_re.fullmatch("c.[1A>G]"))
        self.assertIsNone(multi_variant_re.fullmatch("r.1a>u"))

    def test_does_not_match_mixed_multi(self):
        self.assertIsNone(multi_variant_re.fullmatch("c.[1A>G;Lys4Gly]"))


class TestValidateSingle(TestCase):
    def test_does_not_match_invalid_prefix(self):
        self.assertIsNone(single_variant_re.fullmatch("f.1A>G"))

    def test_can_validate_protein(self):
        self.assertIsNotNone(single_variant_re.fullmatch("p.His4_Gln5insAla"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(His4_Gln5insAla)"))

    def test_can_validate_dna(self):
        self.assertIsNotNone(single_variant_re.fullmatch("c.123A>G"))
        self.assertIsNotNone(single_variant_re.fullmatch("g.19del"))
        self.assertIsNotNone(single_variant_re.fullmatch("m.19_21ins(5)"))
        self.assertIsNotNone(single_variant_re.fullmatch("m.19_21insXXX"))
        self.assertIsNotNone(single_variant_re.fullmatch("n.123_127delinsAAA"))

    def test_can_validate_rna(self):
        self.assertIsNotNone(single_variant_re.fullmatch("r.123a>g"))
        self.assertIsNotNone(single_variant_re.fullmatch("r.19del"))

    def test_does_not_match_invalid_multi(self):
        self.assertIsNone(single_variant_re.fullmatch("c.[1A>G]"))

    def test_does_not_match_invalid(self):
        self.assertIsNone(single_variant_re.fullmatch("c."))


class TestAnyVariantPattern(TestCase):
    def test_matches_single_variant(self):
        self.assertIsNotNone(any_variant_re.fullmatch("c.123A>G"))
        self.assertIsNotNone(any_variant_re.fullmatch("g.19del"))
        self.assertIsNotNone(any_variant_re.fullmatch("m.19_21ins(5)"))
        self.assertIsNotNone(any_variant_re.fullmatch("m.19_21insXXX"))
        self.assertIsNotNone(any_variant_re.fullmatch("n.123_127delinsAAA"))
        self.assertIsNotNone(any_variant_re.fullmatch("r.19del"))
        self.assertIsNotNone(any_variant_re.fullmatch("p.Arg78_Gly79ins23"))
        self.assertIsNotNone(any_variant_re.fullmatch("p.(Arg78_Gly79ins23)"))

    def test_matches_multi_variant(self):
        self.assertIsNotNone(any_variant_re.fullmatch("c.[123A>G;19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("g.[123A>G;19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("m.[123A>G;19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("n.[123A>G;19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("r.[123a>g;19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("r.[123a>g,19del]"))
        self.assertIsNotNone(any_variant_re.fullmatch("p.[Cys28fs;Cys28delinsGly]"))

    def test_none_if_no_match(self):
        self.assertIsNone(any_variant_re.fullmatch("c.123A>F"))
        self.assertIsNone(
            any_variant_re.fullmatch("p.[His4_Gln5insAla,Cys28fs,Cys28delinsVal]")
        )
        self.assertIsNone(any_variant_re.fullmatch("p.[C28fs;(C28delinsG)]"))
        self.assertIsNone(any_variant_re.fullmatch("p.[Cys28fs;C28delinsG]"))
