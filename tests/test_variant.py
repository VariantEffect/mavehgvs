import unittest

from hypothesis import given, strategies as st
from fqfa.constants import AA_CODES

from mavehgvs.exceptions import MaveHgvsParseError
from mavehgvs.variant import Variant
from mavehgvs.position import VariantPosition
from .strategies import (
    AMINO_ACIDS,
    plain_positions,
    utr_intron_positions,
    intron_offset_positions,
    amino_acid_position_strings,
    amino_acid_sequences,
    sequences,
)

# --- strategies for building strings that Variant can fully parse: unlike the
# patterns-level strategies in tests/test_patterns/__init__.py, these respect
# Variant's extra semantic rules (start/end ordering, adjacency for insertions,
# distinct/non-overlapping multi-variant positions).

NUCLEOTIDE_PREFIXES = "cngmor"
ALL_PREFIXES = NUCLEOTIDE_PREFIXES + "p"


def _nucleotide_position_strategy(prefix: str) -> st.SearchStrategy:
    """The position strategy matching the grammar for the given nucleotide
    prefix (c, n, g, m, o, or r)."""
    if prefix == "c":
        return utr_intron_positions()
    elif prefix in "nr":
        return intron_offset_positions()
    else:
        return plain_positions()


def _nucleotide_alphabet(prefix: str) -> str:
    return "acgu" if prefix == "r" else "ACGT"


@st.composite
def sub_variant_strings(draw) -> tuple:
    """A (prefix, string-without-prefix) pair for a substitution variant, valid
    for any of the seven MAVE-HGVS prefixes."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    if prefix == "p":
        pos = draw(amino_acid_position_strings())
        new = draw(st.sampled_from(AMINO_ACIDS))
        return prefix, f"{pos}{new}"
    else:
        pos = draw(_nucleotide_position_strategy(prefix))
        alphabet = _nucleotide_alphabet(prefix)
        ref = draw(st.sampled_from(alphabet))
        new = draw(st.sampled_from(alphabet))
        return prefix, f"{pos}{ref}>{new}"


@st.composite
def fs_strings(draw) -> tuple:
    """A (prefix, string-without-prefix) pair for a protein frameshift variant."""
    return "p", f"{draw(amino_acid_position_strings())}fs"


@st.composite
def single_position_variant_strings(draw, kinds=("del", "dup")) -> tuple:
    """A (prefix, string-without-prefix, kind) triple for a del/dup event using a
    single position, so no start/end ordering is required."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    if prefix == "p":
        pos = draw(amino_acid_position_strings())
    else:
        pos = draw(_nucleotide_position_strategy(prefix))
    kind = draw(st.sampled_from(kinds))
    return prefix, f"{pos}{kind}", kind


@st.composite
def ordered_plain_position_pair(draw) -> tuple:
    """Two distinct plain (non-extended) integer positions, in ascending order."""
    a = draw(st.integers(min_value=1, max_value=10**6))
    b = a + draw(st.integers(min_value=1, max_value=1000))
    return a, b


def _format_ranged_body(start: str, end: str, kind: str, seq) -> str:
    if kind == "delins":
        return f"{start}_{end}{kind}{seq}"
    return f"{start}_{end}{kind}"


@st.composite
def ranged_variant_parts(draw, kinds=("del", "dup", "delins")) -> tuple:
    """A (prefix, start, end, kind, seq) tuple with start/end in ascending
    order (using plain, non-extended positions so ordering is unambiguous).
    ``seq`` is None unless ``kind`` is "delins"."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    a, b = draw(ordered_plain_position_pair())
    if prefix == "p":
        start = f"{draw(st.sampled_from(AMINO_ACIDS))}{a}"
        end = f"{draw(st.sampled_from(AMINO_ACIDS))}{b}"
    else:
        start, end = str(a), str(b)
    kind = draw(st.sampled_from(kinds))
    seq = None
    if kind == "delins":
        if prefix == "p":
            seq = draw(amino_acid_sequences())
        else:
            seq = draw(sequences(_nucleotide_alphabet(prefix)))
    return prefix, start, end, kind, seq


@st.composite
def ranged_variant_strings(draw, kinds=("del", "dup", "delins")) -> tuple:
    """A (prefix, string-without-prefix, kind) triple for a del/dup/delins event
    using an ordered start/end position range."""
    prefix, start, end, kind, seq = draw(ranged_variant_parts(kinds=kinds))
    return prefix, _format_ranged_body(start, end, kind, seq), kind


@st.composite
def adjacent_plain_position_pair(draw) -> tuple:
    """Two adjacent plain integer positions, e.g. (8, 9)."""
    a = draw(st.integers(min_value=1, max_value=10**6))
    return a, a + 1


@st.composite
def non_adjacent_plain_position_pair(draw) -> tuple:
    """Two ordered plain integer positions that are not adjacent (gap >= 2)."""
    a = draw(st.integers(min_value=1, max_value=10**6))
    b = a + draw(st.integers(min_value=2, max_value=1000))
    return a, b


def _ins_string(draw, prefix: str, a: int, b: int) -> str:
    """Build an insertion event string (without the prefix) for the given
    position pair, drawing an appropriate inserted sequence for the prefix."""
    if prefix == "p":
        start = f"{draw(st.sampled_from(AMINO_ACIDS))}{a}"
        end = f"{draw(st.sampled_from(AMINO_ACIDS))}{b}"
        seq = draw(amino_acid_sequences())
    else:
        start, end = str(a), str(b)
        seq = draw(sequences(_nucleotide_alphabet(prefix)))
    return f"{start}_{end}ins{seq}"


@st.composite
def ins_variant_strings(draw) -> tuple:
    """A (prefix, string-without-prefix) pair for an insertion event using an
    adjacent start/end position pair."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    a, b = draw(adjacent_plain_position_pair())
    return prefix, _ins_string(draw, prefix, a, b)


@st.composite
def non_adjacent_ins_strings(draw) -> tuple:
    """A (prefix, string-without-prefix) pair for an insertion event using a
    non-adjacent start/end position pair, which Variant should reject."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    a, b = draw(non_adjacent_plain_position_pair())
    return prefix, _ins_string(draw, prefix, a, b)


@st.composite
def distinct_ascending_positions(draw, min_count: int = 2, max_count: int = 4):
    """A list of distinct plain integer positions in ascending order."""
    n = draw(st.integers(min_value=min_count, max_value=max_count))
    start = draw(st.integers(min_value=1, max_value=10**5))
    gaps = draw(
        st.lists(
            st.integers(min_value=1, max_value=100), min_size=n - 1, max_size=n - 1
        )
    )
    positions = [start]
    for gap in gaps:
        positions.append(positions[-1] + gap)
    return positions


@st.composite
def multi_sub_variant_strings(draw) -> tuple:
    """A (prefix, full_string_with_prefix, count) triple for a multi-variant
    made of substitutions at distinct, ascending plain positions."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    positions = draw(distinct_ascending_positions())
    parts = []
    for pos in positions:
        if prefix == "p":
            aa = draw(st.sampled_from(AMINO_ACIDS))
            new = draw(st.sampled_from(AMINO_ACIDS))
            parts.append(f"{aa}{pos}{new}")
        else:
            alphabet = _nucleotide_alphabet(prefix)
            ref = draw(st.sampled_from(alphabet))
            new = draw(st.sampled_from(alphabet))
            parts.append(f"{pos}{ref}>{new}")
    s = f"{prefix}.[{';'.join(parts)}]"
    return prefix, s, len(positions)


@st.composite
def overlapping_multi_sub_strings(draw) -> str:
    """A multi-variant string with two substitutions at the identical
    position (which Variant should reject as an overlap)."""
    prefix = draw(st.sampled_from(ALL_PREFIXES))
    pos = draw(st.integers(min_value=1, max_value=10**6))
    if prefix == "p":
        aa = draw(st.sampled_from(AMINO_ACIDS))
        new1 = draw(st.sampled_from(AMINO_ACIDS))
        new2 = draw(st.sampled_from(AMINO_ACIDS))
        parts = [f"{aa}{pos}{new1}", f"{aa}{pos}{new2}"]
    else:
        alphabet = _nucleotide_alphabet(prefix)
        ref = draw(st.sampled_from(alphabet))
        new1 = draw(st.sampled_from(alphabet))
        new2 = draw(st.sampled_from(alphabet))
        parts = [f"{pos}{ref}>{new1}", f"{pos}{ref}>{new2}"]
    return f"{prefix}.[{parts[0]};{parts[1]}]"


@st.composite
def dna_target_and_matching_sub(draw) -> tuple:
    """A (target, position, ref, new) tuple where ``ref`` is the actual base at
    ``position`` (1-based) in ``target``."""
    target = draw(st.text(alphabet="ACGT", min_size=1, max_size=30))
    idx = draw(st.integers(min_value=1, max_value=len(target)))
    ref = target[idx - 1]
    new = draw(st.sampled_from([b for b in "ACGT" if b != ref]))
    return target, idx, ref, new


@st.composite
def dna_target_and_mismatching_sub(draw) -> tuple:
    """A (target, position, wrong_ref, new) tuple where ``wrong_ref`` is
    guaranteed to differ from the actual base at ``position`` in ``target``."""
    target = draw(st.text(alphabet="ACGT", min_size=1, max_size=30))
    idx = draw(st.integers(min_value=1, max_value=len(target)))
    actual_ref = target[idx - 1]
    wrong_ref = draw(st.sampled_from([b for b in "ACGT" if b != actual_ref]))
    new = draw(st.sampled_from("ACGT"))
    return target, idx, wrong_ref, new


_PROTEIN_TARGET_LETTERS = [c for c in AA_CODES if c != "*"]


@st.composite
def protein_target_and_matching_sub(draw) -> tuple:
    """A (target, position, aa3, new) tuple where ``aa3`` is the three-letter
    code for the actual residue at ``position`` (1-based) in ``target``."""
    target = "".join(
        draw(
            st.lists(st.sampled_from(_PROTEIN_TARGET_LETTERS), min_size=1, max_size=20)
        )
    )
    idx = draw(st.integers(min_value=1, max_value=len(target)))
    aa3 = AA_CODES[target[idx - 1]]
    new = draw(st.sampled_from(AMINO_ACIDS))
    return target, idx, aa3, new


@st.composite
def protein_target_and_mismatching_sub(draw) -> tuple:
    """A (target, position, wrong_aa3, new) tuple where ``wrong_aa3`` is
    guaranteed to differ from the actual residue at ``position`` in ``target``."""
    target = "".join(
        draw(
            st.lists(st.sampled_from(_PROTEIN_TARGET_LETTERS), min_size=1, max_size=20)
        )
    )
    idx = draw(st.integers(min_value=1, max_value=len(target)))
    actual_letter = target[idx - 1]
    wrong_letter = draw(
        st.sampled_from([c for c in _PROTEIN_TARGET_LETTERS if c != actual_letter])
    )
    wrong_aa3 = AA_CODES[wrong_letter]
    new = draw(st.sampled_from(AMINO_ACIDS))
    return target, idx, wrong_aa3, new


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
            "NM_001301.4:c.122-6T>A",
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
            "NM_001301.4::c.122-6T>A",
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
            "p.Glu12_Gly14=",
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
            *[f"{prefix}.=" for prefix in tuple("gmocnr")],
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
            "NM_002002.3:c.[1_35del;78+5_78+10del;122T>A]",
        ]

        invalid_variant_strings = [
            "p.[Glu27Trp;=;Ter345Lys]",
            "p.[(=);Gly18del;Glu27Trp;Ter345Lys]",
            "p.[Gln7_Asn19=;Glu27Trp;Ter345Lys]",
            "c.[12T>A;=;78+5_78+10del]",
            "c.[1_3=;12T>A;78+5_78+10del]",
            "p.[Glu27fs;Arg48Lys]",
            "p.[Glu27fs;Arg48fs]",
            "NM_002002.3::c.[1_35del;78+5_78+10del;122T>A]",
            "NM_002002.3:c.1_35del;78+5_78+10del;122T>A",
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
                },
                "p.=",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "p",
                    "synonymous": True,
                },
                "p.(=)",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "c",
                },
                "c.=",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "p",
                    "start_position": "27",
                    "start_target": "Glu",
                    "end_position": "27",
                    "end_target": "Glu",
                },
                "p.Glu27=",
            ),
            (
                {
                    "variant_type": "equal",
                    "prefix": "p",
                    "start_position": "12",
                    "start_target": "Glu",
                    "end_position": "14",
                    "end_target": "Gly",
                },
                "p.Glu12_Gly14=",
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
            (
                {
                    "variant_type": "del",
                    "prefix": "p",
                    "start_position": 33,
                    "start_target": "Arg",
                    "end_position": 33,
                    "end_target": "Arg",
                },
                "p.Arg33del",
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
                    "prefix": "c",
                    "start_position": "45",
                    "end_position": "45",
                    "variant": "AGA",
                },
                "c.45delinsAGA",
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
            (
                {
                    "variant_type": "delins",
                    "prefix": "p",
                    "start_position": 50,
                    "start_target": "Arg",
                    "end_position": 50,
                    "end_target": "Arg",
                    "variant": "AlaGly",
                },
                "p.Arg50delinsAlaGly",
            ),
        ]

        invalid_dicts = [
            {
                "variant_type": "delins",
                "prefix": "p",
                "start_position": 50,
                "start_target": "Arg",
                "end_position": 50,
                "end_target": "Cys",
                "variant": "AlaGly",
            },
            {
                "variant_type": "equal",
                "prefix": "p",
                "start_position": "27",
                "start_target": "Glu",
                "end_position": "27",
                "end_target": "Asp",
            },
        ]

        for d, s in valid_dict_tuples:
            with self.subTest(d=d, s=s):
                self.assertEqual(Variant(s), Variant(d))

        for d in invalid_dicts:
            with self.subTest(d=d):
                with self.assertRaises(MaveHgvsParseError):
                    Variant(d)

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
            {
                "variant_type": "fs",
                "prefix": "p",
                "position": 80,
                "target": "Cys",
                "start_position": 23,
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

    def test_components(self):
        variant_strings = [
            ("p.[Glu27Trp;Ter345Lys]", ("p.Glu27Trp", "p.Ter345Lys")),
            ("p.[Glu27Trp;Lys212fs]", ("p.Glu27Trp", "p.Lys212fs")),
            (
                "p.[Gly18del;Glu27Trp;Ter345Lys]",
                ("p.Gly18del", "p.Glu27Trp", "p.Ter345Lys"),
            ),
            (
                "p.[Gln7_Asn19del;Glu27Trp;Ter345Lys]",
                ("p.Gln7_Asn19del", "p.Glu27Trp", "p.Ter345Lys"),
            ),
            (
                "c.[1_35del;78+5_78+10del;122T>A]",
                ("c.1_35del", "c.78+5_78+10del", "c.122T>A"),
            ),
            ("p.Glu27Trp", ("p.Glu27Trp",)),
            ("NP_002002.3:p.Glu27Trp", ("NP_002002.3:p.Glu27Trp",)),
            (
                "NP_002002.3:p.[Glu27Trp;Lys212fs]",
                ("NP_002002.3:p.Glu27Trp", "NP_002002.3:p.Lys212fs"),
            ),
        ]

        for s, expected_components in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(all([c in expected_components for c in v.components()]))


# TODO: multi-variant test cases
class TestMiscProperties(unittest.TestCase):
    def test_prefix(self):
        variant_tuples = [(prefix, f"{prefix}.=") for prefix in tuple("gmocnr")]

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
            ("NM_007294.4", "NM_007294.4:c.[122-6T>A;153C>T]"),
        ]

        for t, s in variant_tuples:
            with self.subTest(t=t, s=s):
                v = Variant(s)
                self.assertEqual(t, v.target_id)

        for _, s in variant_tuples:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))


class TestVariantHypothesisRoundTrip(unittest.TestCase):
    """Property-based tests generalizing the fixed round-trip examples above
    across all seven MAVE-HGVS prefixes."""

    @given(pv=sub_variant_strings())
    def test_sub_round_trip(self, pv: tuple) -> None:
        prefix, body = pv
        s = f"{prefix}.{body}"
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertEqual("sub", v.variant_type)
        self.assertEqual(prefix, v.prefix)

    @given(pv=fs_strings())
    def test_fs_round_trip(self, pv: tuple) -> None:
        prefix, body = pv
        s = f"{prefix}.{body}"
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertEqual("fs", v.variant_type)

    @given(pvk=single_position_variant_strings())
    def test_single_position_del_dup_round_trip(self, pvk: tuple) -> None:
        prefix, body, kind = pvk
        s = f"{prefix}.{body}"
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertEqual(kind, v.variant_type)

    @given(pvk=ranged_variant_strings())
    def test_ranged_del_dup_delins_round_trip(self, pvk: tuple) -> None:
        prefix, body, kind = pvk
        s = f"{prefix}.{body}"
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertEqual(kind, v.variant_type)

    @given(pv=ins_variant_strings())
    def test_ins_round_trip(self, pv: tuple) -> None:
        prefix, body = pv
        s = f"{prefix}.{body}"
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertEqual("ins", v.variant_type)

    @given(pv=non_adjacent_ins_strings())
    def test_non_adjacent_ins_rejected(self, pv: tuple) -> None:
        prefix, body = pv
        with self.assertRaises(MaveHgvsParseError):
            Variant(f"{prefix}.{body}")


class TestVariantHypothesisOrdering(unittest.TestCase):
    @given(parts=ranged_variant_parts(kinds=("del", "dup")))
    def test_relaxed_ordering_swaps_positions(self, parts: tuple) -> None:
        prefix, start, end, kind, seq = parts
        canonical_s = f"{prefix}.{_format_ranged_body(start, end, kind, seq)}"
        reversed_s = f"{prefix}.{_format_ranged_body(end, start, kind, seq)}"

        with self.assertRaises(MaveHgvsParseError):
            Variant(reversed_s)

        v = Variant(reversed_s, relaxed_ordering=True)
        self.assertEqual(canonical_s, str(v))


class TestVariantHypothesisMultiVariant(unittest.TestCase):
    @given(pv=multi_sub_variant_strings())
    def test_multi_variant_round_trip(self, pv: tuple) -> None:
        prefix, s, count = pv
        v = Variant(s)
        self.assertEqual(s, str(v))
        self.assertTrue(v.is_multi_variant())
        self.assertEqual(count, v.variant_count)

    @given(s=overlapping_multi_sub_strings())
    def test_multi_variant_same_position_rejected(self, s: str) -> None:
        with self.assertRaises(MaveHgvsParseError):
            Variant(s)


class TestVariantHypothesisTargetSequence(unittest.TestCase):
    @given(t=dna_target_and_matching_sub())
    def test_dna_sub_matches_target(self, t: tuple) -> None:
        target, idx, ref, new = t
        s = f"c.{idx}{ref}>{new}"
        v = Variant(s, targetseq=target)
        self.assertEqual(s, str(v))

    @given(t=dna_target_and_mismatching_sub())
    def test_dna_sub_mismatched_ref_rejected(self, t: tuple) -> None:
        target, idx, wrong_ref, new = t
        s = f"c.{idx}{wrong_ref}>{new}"
        with self.assertRaises(MaveHgvsParseError):
            Variant(s, targetseq=target)

    @given(
        target=st.text(alphabet="ACGT", min_size=1, max_size=30),
        extra=st.integers(min_value=1, max_value=1000),
        ref=st.sampled_from("ACGT"),
        new=st.sampled_from("ACGT"),
    )
    def test_dna_sub_out_of_bounds_rejected(
        self, target: str, extra: int, ref: str, new: str
    ) -> None:
        idx = len(target) + extra
        s = f"c.{idx}{ref}>{new}"
        with self.assertRaises(MaveHgvsParseError):
            Variant(s, targetseq=target)

    @given(t=protein_target_and_matching_sub())
    def test_protein_sub_matches_target(self, t: tuple) -> None:
        target, idx, aa3, new = t
        s = f"p.{aa3}{idx}{new}"
        v = Variant(s, targetseq=target)
        self.assertEqual(s, str(v))

    @given(t=protein_target_and_mismatching_sub())
    def test_protein_sub_mismatched_target_rejected(self, t: tuple) -> None:
        target, idx, wrong_aa3, new = t
        s = f"p.{wrong_aa3}{idx}{new}"
        with self.assertRaises(MaveHgvsParseError):
            Variant(s, targetseq=target)


if __name__ == "__main__":
    unittest.main()
