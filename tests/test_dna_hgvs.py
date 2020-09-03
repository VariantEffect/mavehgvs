from unittest import TestCase

from hgvsp.dna import (
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
        self.assertIsNotNone(substitution_re.fullmatch("123A>G"))
        for c in ("c", "n", "g", "m"):
            self.assertIsNotNone(substitution_re.fullmatch(f"{c}.123A>G"))
        self.assertIsNotNone(substitution_re.fullmatch("123A>C"))
        self.assertIsNotNone(substitution_re.fullmatch("123A>N"))
        self.assertIsNotNone(substitution_re.fullmatch("*123A>G"))
        self.assertIsNotNone(substitution_re.fullmatch("-123A>G"))
        self.assertIsNotNone(substitution_re.fullmatch("-123+45A>G"))
        self.assertIsNotNone(substitution_re.fullmatch("*123-45A>G"))
        self.assertIsNotNone(substitution_re.fullmatch("93+1G>T"))
        self.assertIsNotNone(substitution_re.fullmatch("54G>H"))
        self.assertIsNotNone(substitution_re.fullmatch("54="))
        self.assertIsNotNone(substitution_re.fullmatch("54=/T>C"))
        self.assertIsNotNone(substitution_re.fullmatch("54=//T>C"))

    def test_does_not_match_invalid_substitutions(self):
        self.assertIsNone(substitution_re.fullmatch(""))
        self.assertIsNone(substitution_re.fullmatch("A>G"))
        self.assertIsNone(substitution_re.fullmatch("*A>G"))
        self.assertIsNone(substitution_re.fullmatch("12A=G"))
        self.assertIsNone(substitution_re.fullmatch("12A>E"))
        self.assertIsNone(substitution_re.fullmatch("12A<E"))
        self.assertIsNone(substitution_re.fullmatch("12>A"))
        self.assertIsNone(substitution_re.fullmatch("+12A>G"))
        self.assertIsNone(substitution_re.fullmatch("123A>X"))

    def test_valid_deletions_pass(self):
        self.assertIsNotNone(deletion_re.fullmatch("19del"))
        for c in ("c", "n", "g", "m"):
            self.assertIsNotNone(deletion_re.fullmatch(f"{c}.123delA"))
        self.assertIsNotNone(deletion_re.fullmatch("19delT"))
        self.assertIsNotNone(deletion_re.fullmatch("19_21del"))
        self.assertIsNotNone(deletion_re.fullmatch("*183_186+48del"))
        self.assertIsNotNone(deletion_re.fullmatch("*183+45_186del"))
        self.assertIsNotNone(deletion_re.fullmatch("1704+1del"))
        self.assertIsNotNone(deletion_re.fullmatch("4072-1234_5155-246del"))
        self.assertIsNotNone(
            deletion_re.fullmatch("(4071+1_4072-1)_(5154+1_5155-1)del")
        )
        self.assertIsNotNone(deletion_re.fullmatch("720_991del"))
        self.assertIsNotNone(deletion_re.fullmatch("(?_-245)_(31+1_32-1)del"))
        self.assertIsNotNone(deletion_re.fullmatch("(?_-1)_(*1_?)del"))
        self.assertIsNotNone(deletion_re.fullmatch("19_21=/del"))
        self.assertIsNotNone(deletion_re.fullmatch("19_21del=//del"))

    def test_does_not_match_invalid_deletions(self):
        self.assertIsNone(deletion_re.fullmatch("19delE"))
        self.assertIsNone(deletion_re.fullmatch(""))
        self.assertIsNone(deletion_re.fullmatch("delA"))
        self.assertIsNone(deletion_re.fullmatch("4071+1_4072-1_5154+1_5155-1del"))
        self.assertIsNone(deletion_re.fullmatch("(?_-1)_(+1_?)del"))
        self.assertIsNone(deletion_re.fullmatch("1704+1delAAA"))
        self.assertIsNone(deletion_re.fullmatch("19_21del(5)"))
        self.assertIsNone(deletion_re.fullmatch("19_21delTTT"))

    def test_valid_insertions_pass(self):
        self.assertIsNotNone(insertion_re.fullmatch("169_170insA"))
        self.assertIsNotNone(insertion_re.fullmatch("240_241insAGG"))
        self.assertIsNotNone(insertion_re.fullmatch("761_762insNNNNN"))
        self.assertIsNotNone(insertion_re.fullmatch("32717298_32717299ins(100)"))
        self.assertIsNotNone(insertion_re.fullmatch("761_762insN"))
        self.assertIsNotNone(insertion_re.fullmatch("(222_226)insG"))
        for c in ("c", "n", "g", "m"):
            self.assertIsNotNone(insertion_re.fullmatch(f"{c}.123_124insA"))

    def test_does_not_match_invalid_insertions(self):
        self.assertIsNone(insertion_re.fullmatch("19insR"))
        self.assertIsNone(insertion_re.fullmatch(""))
        self.assertIsNone(insertion_re.fullmatch("insA"))
        self.assertIsNone(insertion_re.fullmatch("(4071+1_4072)-(1_5154+1_5155-1)ins"))
        self.assertIsNone(insertion_re.fullmatch("(?_-1)_(+1_?)ins"))
        self.assertIsNone(insertion_re.fullmatch("1704+1insAAA"))

    def test_valid_delins_passes(self):
        self.assertIsNotNone(delins_re.fullmatch("6775delinsGA"))
        self.assertIsNotNone(delins_re.fullmatch("6775_6777delinsC"))
        self.assertIsNotNone(delins_re.fullmatch("?_6777delinsC"))
        self.assertIsNotNone(delins_re.fullmatch("*?_45+1delinsC"))
        self.assertIsNotNone(delins_re.fullmatch("?_?delinsC"))
        self.assertIsNotNone(delins_re.fullmatch("142_144delinsTGG"))
        self.assertIsNotNone(delins_re.fullmatch("9002_9009delinsTTT"))
        self.assertIsNotNone(delins_re.fullmatch("9002_9009delins(5)"))
        for c in ("c", "n", "g", "m"):
            self.assertIsNotNone(delins_re.fullmatch(f"{c}.123_127delinsA"))

    def test_does_not_match_invalid_delins(self):
        self.assertIsNone(delins_re.fullmatch("19delinsE"))
        self.assertIsNone(delins_re.fullmatch(""))
        self.assertIsNone(delins_re.fullmatch("delinsA"))
        self.assertIsNone(delins_re.fullmatch("(4071+1_4072)-(1_5154+1_5155-1)delins"))
        self.assertIsNone(delins_re.fullmatch("(?_-1)_(+1_?)delins"))

    def test_any_event_re_matches_any(self):
        self.assertIsNotNone(any_event_re.fullmatch("6775delinsGA"))
        self.assertIsNotNone(any_event_re.fullmatch("169_170insA"))
        self.assertIsNotNone(any_event_re.fullmatch("19del"))
        self.assertIsNotNone(any_event_re.fullmatch("*123-45A>G"))


class TestVariantRegexPatterns(TestCase):
    def test_single_var_re_matches_each_variant_type(self):
        self.assertIsNotNone(single_variant_re.fullmatch("c.123A>G"))
        self.assertIsNotNone(
            single_variant_re.fullmatch("g.(4071+1_4072-1)_(5154+1_5155-1)del")
        )
        self.assertIsNotNone(single_variant_re.fullmatch("n.240_241insAGG"))
        self.assertIsNotNone(single_variant_re.fullmatch("m.9002_9009delins(5)"))

    def test_multi_var_re_matches_multi_variants(self):
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "c.[123A>G;19del;" "240_241insAGG;9002_9009delins(5)]"
            )
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "g.[123=/A>G;19delT;" "240_241insAGG;9002_9009delins(5)]"
            )
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "m.[123=;(?_-1)_(*1_?)del;"
                "(?_-245)_(31+1_32-1)del;9002_9009delinsGGG]"
            )
        )
        self.assertIsNotNone(
            multi_variant_re.fullmatch(
                "n.[54=//T>C;*183+45_186del;" "32717298_32717299ins(100);6775delinsGA]"
            )
        )

        # Non-multi should be none
        self.assertIsNone(multi_variant_re.fullmatch("c.[123=;]"))
        self.assertIsNone(multi_variant_re.fullmatch("c.[123A>G]"))
