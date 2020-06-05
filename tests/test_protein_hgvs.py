from unittest import TestCase

from hgvsp.protein import (
    substitution_re,
    deletion_re,
    delins_re,
    insertion_re,
    frame_shift_re,
    single_variant_re,
    multi_variant_re,
    any_event_re,
)


class TestEventValidators(TestCase):
    def test_valid_substitutions_pass(self):
        self.assertIsNotNone(substitution_re.fullmatch("Trp24Cys"))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24???"))
        self.assertIsNotNone(substitution_re.fullmatch("Cys188="))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24*"))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24Ter"))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24Ter^Ala^G"))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24?"))
        self.assertIsNotNone(substitution_re.fullmatch("Trp24=/Cys"))
        self.assertIsNotNone(substitution_re.fullmatch("p.Trp24=/Cys"))
        self.assertIsNotNone(substitution_re.fullmatch("p.(Trp24=/Cys)"))
        self.assertIsNotNone(substitution_re.fullmatch("0"))
        self.assertIsNotNone(substitution_re.fullmatch("?"))
        self.assertIsNotNone(substitution_re.fullmatch("p.="))

    def test_does_not_match_invalid_substitutions(self):
        self.assertIsNone(substitution_re.fullmatch(""))
        self.assertIsNone(substitution_re.fullmatch("a>g"))
        self.assertIsNone(substitution_re.fullmatch("*a>g"))
        self.assertIsNone(substitution_re.fullmatch("1a>a"))
        self.assertIsNone(substitution_re.fullmatch("12a=g"))
        self.assertIsNone(substitution_re.fullmatch("12a>E"))
        self.assertIsNone(substitution_re.fullmatch("12a<E"))
        self.assertIsNone(substitution_re.fullmatch("12-1>a"))
        self.assertIsNone(substitution_re.fullmatch("+12a>g"))

    def test_valid_deletions_pass(self):
        self.assertIsNotNone(deletion_re.fullmatch("Val7del"))
        self.assertIsNotNone(deletion_re.fullmatch("Lys23_Val25del"))
        self.assertIsNotNone(deletion_re.fullmatch("Trp4del"))
        self.assertIsNotNone(deletion_re.fullmatch("???4del"))
        self.assertIsNotNone(deletion_re.fullmatch("?4del"))
        self.assertIsNotNone(deletion_re.fullmatch("Gly2_Met46del"))
        self.assertIsNotNone(deletion_re.fullmatch("Val7=/del"))
        self.assertIsNotNone(deletion_re.fullmatch("p.Val7=/del"))
        self.assertIsNotNone(deletion_re.fullmatch("p.(Val7=/del)"))

    def test_does_not_match_invalid_deletions(self):
        self.assertIsNone(deletion_re.fullmatch("19delR"))
        self.assertIsNone(deletion_re.fullmatch(""))
        self.assertIsNone(deletion_re.fullmatch("dela"))
        self.assertIsNone(deletion_re.fullmatch("4071+1_4072-1_5154+1_5155-1del"))
        self.assertIsNone(deletion_re.fullmatch("(?_-1)_(+1_?)del"))
        self.assertIsNone(deletion_re.fullmatch("1704+1delaaa"))
        self.assertIsNone(deletion_re.fullmatch("19_21del(5)"))
        self.assertIsNone(deletion_re.fullmatch("19_21deluuu"))

    def test_valid_insertions_pass(self):
        self.assertIsNotNone(insertion_re.fullmatch("His4_Gln5insAla"))
        self.assertIsNotNone(insertion_re.fullmatch("His4_Gln5ins???"))
        self.assertIsNotNone(insertion_re.fullmatch("His4_Gln5ins?"))
        self.assertIsNotNone(insertion_re.fullmatch("His4_Gln5insAla^Gly^Ser"))
        self.assertIsNotNone(insertion_re.fullmatch("Lys2_Gly3insGlnSerLys"))
        self.assertIsNotNone(insertion_re.fullmatch("Met3_His4insGlyTer"))
        self.assertIsNotNone(insertion_re.fullmatch("Arg78_Gly79ins23"))
        self.assertIsNotNone(insertion_re.fullmatch("Ser332_Ser333ins(1)"))
        self.assertIsNotNone(insertion_re.fullmatch("Val582_Asn583ins(5)"))
        self.assertIsNotNone(insertion_re.fullmatch("Val582_Asn583insX"))
        self.assertIsNotNone(insertion_re.fullmatch("Val582_Asn583insXXXXX"))
        self.assertIsNotNone(insertion_re.fullmatch("p.Val582_Asn583ins(5)"))
        self.assertIsNotNone(insertion_re.fullmatch("p.(Val582_Asn583ins(5))"))

    def test_does_not_match_invalid_insertions(self):
        self.assertIsNone(insertion_re.fullmatch("Val582insXXXXX"))
        self.assertIsNone(insertion_re.fullmatch(""))
        self.assertIsNone(insertion_re.fullmatch("ins"))
        self.assertIsNone(insertion_re.fullmatch("(Val582_Asn583)insXXXXX"))
        self.assertIsNone(insertion_re.fullmatch("Val582_Asn583ins"))
        self.assertIsNone(insertion_re.fullmatch("Val582+1_Asn583insXXXXX"))

    def test_valid_delins_passes(self):
        self.assertIsNotNone(delins_re.fullmatch("Cys28delinsTrpVal"))
        self.assertIsNotNone(delins_re.fullmatch("C28_L29delinsT"))
        self.assertIsNotNone(delins_re.fullmatch("C28_L29delins*"))
        self.assertIsNotNone(delins_re.fullmatch("Cys28delinsTrpVal"))
        self.assertIsNotNone(delins_re.fullmatch("Cys28delins???Val"))
        self.assertIsNotNone(delins_re.fullmatch("Cys28delins?Val"))
        self.assertIsNotNone(
            delins_re.fullmatch("Glu125_Ala132delinsGlyLeuHisArgPheIleValLeu")
        )
        self.assertIsNotNone(delins_re.fullmatch("C28_L29delinsT^G^L"))
        self.assertIsNotNone(delins_re.fullmatch("p.C28_L29delinsT^G^L"))
        self.assertIsNotNone(delins_re.fullmatch("p.(C28_L29delinsT^G^L)"))

    def test_does_not_match_invalid_delins(self):
        self.assertIsNone(delins_re.fullmatch("Cys28delinsJ"))
        self.assertIsNone(delins_re.fullmatch(""))
        self.assertIsNone(delins_re.fullmatch("(Cys28_Cys)delinsTrpVal"))
        self.assertIsNone(delins_re.fullmatch("C28_L29delinsTG^G^L"))
        self.assertIsNone(delins_re.fullmatch("Cys28+5delinsZ"))
        self.assertIsNone(delins_re.fullmatch("*?_45+1delinsg"))

    def test_valid_frameshift_passes(self):
        self.assertIsNotNone(frame_shift_re.fullmatch("Arg97ProfsTer23"))
        self.assertIsNotNone(frame_shift_re.fullmatch("Glu5ValfsTer5"))
        self.assertIsNotNone(frame_shift_re.fullmatch("Ile327Argfs*?"))
        self.assertIsNotNone(frame_shift_re.fullmatch("Ile327Argfs???5"))
        self.assertIsNotNone(frame_shift_re.fullmatch("Ile327fs"))
        self.assertIsNotNone(frame_shift_re.fullmatch("Gln151Thrfs*9"))
        self.assertIsNotNone(frame_shift_re.fullmatch("p.Ile327fs"))
        self.assertIsNotNone(frame_shift_re.fullmatch("p.(Ile327fs)"))

    def test_does_not_match_invalid_frameshift(self):
        self.assertIsNone(frame_shift_re.fullmatch("Arg97ProfsTer23Pro"))
        self.assertIsNone(frame_shift_re.fullmatch(""))
        self.assertIsNone(frame_shift_re.fullmatch("fsTer"))
        self.assertIsNone(frame_shift_re.fullmatch("Glu5_Val7fsTer5"))
        self.assertIsNone(frame_shift_re.fullmatch("Ile327Argfs*?Ter"))
        self.assertIsNone(frame_shift_re.fullmatch("Ile327fs(4)"))
        self.assertIsNone(frame_shift_re.fullmatch("*?_45+1delinsc"))

    def test_any_event_re_matches_any(self):
        self.assertIsNotNone(any_event_re.fullmatch("Trp24Ter^Ala^G"))
        self.assertIsNotNone(any_event_re.fullmatch("Arg97ProfsTer23"))
        self.assertIsNotNone(any_event_re.fullmatch("Cys28delinsTrpVal"))
        self.assertIsNotNone(any_event_re.fullmatch("Arg78_Gly79ins23"))
        self.assertIsNotNone(any_event_re.fullmatch("Arg78_???79ins23"))
        self.assertIsNotNone(any_event_re.fullmatch("Arg78_X79ins23"))
        self.assertIsNotNone(any_event_re.fullmatch("Trp24=/Cys"))


class TestVariantRegexPatterns(TestCase):
    def test_single_var_re_matches_each_variant_type(self):
        self.assertIsNotNone(single_variant_re.fullmatch("p.Trp24Cys^Gly"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.Trp24X"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(Trp24Cys)"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.Lys23_Val25del"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(Lys23_Val25del)"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.His4_Gln5insAla"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(His4_Gln5insAla)"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.Cys28delinsVal"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(Cys28delinsVal)"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.Cys28fs"))
        self.assertIsNotNone(single_variant_re.fullmatch("p.(Cys28fs)"))

    def test_multi_var_re_matches_multi_variants(self):
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "p.[Trp24Cys;Lys23_Val25del;His4_Gln5insAla;" "Cys28fs;Cys28delinsVal]"
            )
        )
        self.assertIsNotNone(multi_variant_re.fullmatch("p.[(Trp24Cys);(Trp24Cys)]"))
        self.assertIsNotNone(
            multi_variant_re.fullmatch("p.[(Trp24Cys);Cys28delinsVal]")
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch("p.[(His4_Gln5insAla);(Cys28fs);Cys28delinsVal]")
        )

        # Non-multi should be none
        self.assertIsNone(multi_variant_re.fullmatch("p.[Trp24Cys;]"))
        self.assertIsNone(multi_variant_re.fullmatch("p.[(Trp24Cys)]"))
        self.assertIsNone(multi_variant_re.fullmatch("p.[Trp24Cys,]"))
        self.assertIsNone(multi_variant_re.fullmatch("p.[Trp24Cys]"))
        self.assertIsNone(multi_variant_re.fullmatch("p.[Trp24???]"))
