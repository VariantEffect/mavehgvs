from unittest import TestCase

from hgvsp.rna import (
    deletion_re,
    insertion_re,
    delins_re,
    substitution_re,
    single_variant_re,
    multi_variant_re,
    any_event_re,
)


class TestEventValidators(TestCase):
    def test_valid_substitutions_pass(self):
        self.assertIsNotNone(substitution_re.fullmatch("123a>g"))
        self.assertIsNotNone(substitution_re.fullmatch("r.123a>g"))
        self.assertIsNotNone(substitution_re.fullmatch("123a>x"))
        self.assertIsNotNone(substitution_re.fullmatch("123a>n"))
        self.assertIsNotNone(substitution_re.fullmatch("54g>h"))
        self.assertIsNotNone(substitution_re.fullmatch("54="))
        self.assertIsNotNone(substitution_re.fullmatch("54=/u>c"))
        self.assertIsNotNone(substitution_re.fullmatch("54=//u>c"))
        self.assertIsNotNone(substitution_re.fullmatch("0"))
        self.assertIsNotNone(substitution_re.fullmatch("?"))
        self.assertIsNotNone(substitution_re.fullmatch("spl"))

    def test_does_not_match_invalid_substitutions(self):
        self.assertIsNone(substitution_re.fullmatch(""))
        self.assertIsNone(substitution_re.fullmatch("a>g"))
        self.assertIsNone(substitution_re.fullmatch("*a>g"))
        self.assertIsNone(substitution_re.fullmatch("12a=g"))
        self.assertIsNone(substitution_re.fullmatch("12a>E"))
        self.assertIsNone(substitution_re.fullmatch("12a<E"))
        self.assertIsNone(substitution_re.fullmatch("12-1>a"))
        self.assertIsNone(substitution_re.fullmatch("+12a>g"))

    def test_valid_deletions_pass(self):
        self.assertIsNotNone(deletion_re.fullmatch("10del"))
        self.assertIsNotNone(deletion_re.fullmatch("6_8del"))
        self.assertIsNotNone(deletion_re.fullmatch("19_21del"))
        self.assertIsNotNone(deletion_re.fullmatch("(4072_5145)del"))
        self.assertIsNotNone(deletion_re.fullmatch("=/6_8del"))
        self.assertIsNotNone(deletion_re.fullmatch("1704del"))
        self.assertIsNotNone(deletion_re.fullmatch("r.1704del"))

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
        self.assertIsNotNone(insertion_re.fullmatch("426_427insa"))
        self.assertIsNotNone(insertion_re.fullmatch("r.426_427insa"))
        self.assertIsNotNone(insertion_re.fullmatch("756_757insuacu"))
        self.assertIsNotNone(insertion_re.fullmatch("(222_226)insg"))
        self.assertIsNotNone(insertion_re.fullmatch("549_550insn"))
        self.assertIsNotNone(insertion_re.fullmatch("761_762insnnnnn"))
        self.assertIsNotNone(insertion_re.fullmatch("761_762ins(5)"))
        self.assertIsNotNone(
            insertion_re.fullmatch("2949_2950ins[2950-30_2950-12;2950-4_2950-1]")
        )

    def test_does_not_match_invalid_insertions(self):
        self.assertIsNone(insertion_re.fullmatch("19insR"))
        self.assertIsNone(insertion_re.fullmatch(""))
        self.assertIsNone(insertion_re.fullmatch("insa"))
        self.assertIsNone(insertion_re.fullmatch("(4071+1_4072)-(1_5154+1_5155-1)ins"))
        self.assertIsNone(insertion_re.fullmatch("(?_-1)_(+1_?)ins"))
        self.assertIsNone(insertion_re.fullmatch("1704+1insaaa"))

    def test_valid_delins_passes(self):
        self.assertIsNotNone(delins_re.fullmatch("6775delinsga"))
        self.assertIsNotNone(delins_re.fullmatch("r.6775delinsga"))
        self.assertIsNotNone(delins_re.fullmatch("6775_6777delinsc"))
        self.assertIsNotNone(delins_re.fullmatch("?_6777delinsc"))
        self.assertIsNotNone(delins_re.fullmatch("?_?delinsc"))
        self.assertIsNotNone(delins_re.fullmatch("142_144delinsugg"))
        self.assertIsNotNone(delins_re.fullmatch("9002_9009delinsuuu"))
        self.assertIsNotNone(delins_re.fullmatch("9002_9009delins(5)"))

    def test_does_not_match_invalid_delins(self):
        self.assertIsNone(delins_re.fullmatch("19delinsR"))
        self.assertIsNone(delins_re.fullmatch(""))
        self.assertIsNone(delins_re.fullmatch("delinsa"))
        self.assertIsNone(delins_re.fullmatch("(4071+1_4072)-(1_5154+1_5155-1)delins"))
        self.assertIsNone(delins_re.fullmatch("(?_-1)_(+1_?)delins"))
        self.assertIsNone(delins_re.fullmatch("*?_45+1delinsc"))

    def test_any_event_re_matches_any(self):
        self.assertIsNotNone(any_event_re.fullmatch("6775delinsga"))
        self.assertIsNotNone(any_event_re.fullmatch("426_427insa"))
        self.assertIsNotNone(any_event_re.fullmatch("19del"))
        self.assertIsNotNone(any_event_re.fullmatch("123a>g"))


class TestVariantRegexPatterns(TestCase):
    def test_single_var_re_matches_each_variant_type(self):
        self.assertIsNotNone(single_variant_re.fullmatch("r.123a>g"))
        self.assertIsNotNone(single_variant_re.fullmatch("r.=/6_8del"))
        self.assertIsNotNone(
            single_variant_re.fullmatch("r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]")
        )
        self.assertIsNotNone(single_variant_re.fullmatch("r.9002_9009delins(5)"))

    def test_multi_var_re_matches_multi_variants(self):
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "r.[123a>g;19del;"
                "2949_2950ins[2950-30_2950-12;2950-4_2950-1];"
                "9002_9009delins(5)]"
            )
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "r.[123a>g,19del,"
                "2949_2950ins[2950-30_2950-12;2950-4_2950-1],"
                "9002_9009delins(5)]"
            )
        )

        # Non-multi should be none
        self.assertIsNone(multi_variant_re.fullmatch("r.[123=;]"))
        self.assertIsNone(multi_variant_re.fullmatch("r.[123=,]"))
        self.assertIsNone(multi_variant_re.fullmatch("r.[123a>g]"))
