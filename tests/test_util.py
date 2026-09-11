import unittest
import string

from hypothesis import given, strategies as st

from mavehgvs.util import parse_variant_strings
from mavehgvs.variant import Variant
from .test_variant import (
    ALL_PREFIXES,
    sub_variant_strings,
    dna_target_and_matching_sub,
    dna_target_and_mismatching_sub,
)


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


@st.composite
def sub_variant_and_wrong_prefix(draw) -> tuple:
    """A (prefix, body, wrong_prefix) triple where wrong_prefix is a valid
    MAVE-HGVS prefix that differs from the variant's actual prefix."""
    prefix, body = draw(sub_variant_strings())
    wrong_prefix = draw(st.sampled_from(ALL_PREFIXES).filter(lambda p: p != prefix))
    return prefix, body, wrong_prefix


class TestParseVariantStringsHypothesis(unittest.TestCase):
    """Property-based tests generalizing the fixed examples above: any
    generated valid variant string should be wrapped by parse_variant_strings
    exactly the way constructing Variant directly would behave."""

    @given(pv=sub_variant_strings())
    def test_valid_items_match_direct_variant_construction(self, pv: tuple) -> None:
        prefix, body = pv
        s = f"{prefix}.{body}"
        valid, invalid = parse_variant_strings([s])
        self.assertEqual(Variant(s), valid[0])
        self.assertIsNone(invalid[0])

    @given(
        s=st.text(
            alphabet=string.ascii_letters + string.digits + ".>*+-_;[]=", max_size=20
        )
    )
    def test_never_raises_unexpected_exceptions(self, s: str) -> None:
        valid, invalid = parse_variant_strings([s])
        self.assertEqual(1, len(valid))
        self.assertEqual(1, len(invalid))
        self.assertEqual(valid[0] is None, invalid[0] is not None)
        if valid[0] is not None:
            self.assertIsInstance(valid[0], Variant)
        else:
            self.assertIsInstance(invalid[0], str)

    @given(pvw=sub_variant_and_wrong_prefix())
    def test_expected_prefix_filtering(self, pvw: tuple) -> None:
        prefix, body, wrong_prefix = pvw
        s = f"{prefix}.{body}"

        valid, invalid = parse_variant_strings([s], expected_prefix=prefix)
        self.assertIsInstance(valid[0], Variant)
        self.assertIsNone(invalid[0])

        valid, invalid = parse_variant_strings([s], expected_prefix=wrong_prefix)
        self.assertIsNone(valid[0])
        self.assertEqual("unexpected variant prefix", invalid[0])

    @given(
        p=st.text(alphabet=string.printable, min_size=1, max_size=1).filter(
            lambda c: c not in "cgmnopr"
        )
    )
    def test_invalid_expected_prefix_raises_valueerror(self, p: str) -> None:
        with self.assertRaises(ValueError):
            parse_variant_strings(["p.Glu27Trp"], expected_prefix=p)

    @given(t=dna_target_and_matching_sub())
    def test_targetseq_matching_sub(self, t: tuple) -> None:
        target, idx, ref, new = t
        s = f"c.{idx}{ref}>{new}"
        valid, invalid = parse_variant_strings([s], targetseq=target)
        self.assertIsInstance(valid[0], Variant)
        self.assertIsNone(invalid[0])

    @given(t=dna_target_and_mismatching_sub())
    def test_targetseq_mismatching_sub(self, t: tuple) -> None:
        target, idx, wrong_ref, new = t
        s = f"c.{idx}{wrong_ref}>{new}"
        valid, invalid = parse_variant_strings([s], targetseq=target)
        self.assertIsNone(valid[0])
        self.assertIsInstance(invalid[0], str)

    @given(pvs=st.lists(sub_variant_strings(), min_size=1, max_size=5))
    def test_batch_length_and_pairing_invariant(self, pvs: list) -> None:
        strings = [f"{prefix}.{body}" for prefix, body in pvs]
        valid, invalid = parse_variant_strings(strings)
        self.assertEqual(len(strings), len(valid))
        self.assertEqual(len(strings), len(invalid))
        for s, v, i in zip(strings, valid, invalid):
            self.assertEqual(v is None, i is not None)
            self.assertEqual(Variant(s), v)
