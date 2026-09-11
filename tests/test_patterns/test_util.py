import re
import string
import unittest

from hypothesis import assume, given, strategies as st

from mavehgvs.patterns.util import combine_patterns, remove_named_groups


@st.composite
def group_names(draw) -> str:
    first = draw(st.sampled_from(string.ascii_letters + "_"))
    rest = draw(
        st.text(alphabet=string.ascii_letters + string.digits + "_", max_size=8)
    )
    return first + rest


@st.composite
def named_group_pattern_and_example(draw) -> tuple:
    """A pattern of the form (?P<name>(?P<inner>[ch])) together with an example
    string it matches."""
    name = draw(group_names())
    inner = draw(group_names())
    ch = draw(st.sampled_from(string.ascii_lowercase))
    pattern = f"(?P<{name}>(?P<{inner}>[{ch}]))"
    return pattern, ch


class TestCombinePatterns(unittest.TestCase):
    def test_without_groupname(self):
        pattern_tuples = [
            (
                ("(?P<a>(?P<x>[1-9]))", "(?P<b>(?P<x>[1-9]))"),
                "(?:(?P<a>(?P<a_x>[1-9]))|(?P<b>(?P<b_x>[1-9])))",
            )
        ]

        for p1, p2 in pattern_tuples:
            with self.subTest(p1=p1, p2=p2):
                self.assertEqual(combine_patterns(p1), p2)

    def test_with_groupname(self):
        pattern_tuples = [
            (
                ("(?P<a>(?P<x>[1-9]))", "(?P<b>(?P<x>[1-9]))"),
                "test",
                "(?P<test>(?P<a>(?P<a_x>[1-9]))|(?P<b>(?P<b_x>[1-9])))",
            )
        ]

        for p1, g, p2 in pattern_tuples:
            with self.subTest(p1=p1, g=g, p2=p2):
                self.assertEqual(combine_patterns(p1, groupname=g), p2)


class TestRemoveNamedGroups(unittest.TestCase):
    def test_noncapturing(self):
        pattern_tuples = [("(?P<a>(?P<x>[1-9]))", "(?:(?:[1-9]))")]

        for p1, p2 in pattern_tuples:
            with self.subTest(p1=p1, p2=p2):
                self.assertEqual(remove_named_groups(p1, noncapturing=True), p2)

    def test_capturing(self):
        pattern_tuples = [("(?P<a>(?P<x>[1-9]))", "(([1-9]))")]

        for p1, p2 in pattern_tuples:
            with self.subTest(p1=p1, p2=p2):
                self.assertEqual(remove_named_groups(p1, noncapturing=False), p2)


class TestHypothesisRemoveNamedGroups(unittest.TestCase):
    @given(pg=named_group_pattern_and_example(), noncapturing=st.booleans())
    def test_no_named_groups_remain(self, pg: tuple, noncapturing: bool) -> None:
        pattern, example = pg
        stripped = remove_named_groups(pattern, noncapturing=noncapturing)
        self.assertNotIn("(?P<", stripped)
        self.assertIsNotNone(re.compile(stripped).fullmatch(example))


class TestHypothesisCombinePatterns(unittest.TestCase):
    @given(
        name1=group_names(),
        name2=group_names(),
        ch1=st.sampled_from(string.ascii_lowercase),
        ch2=st.sampled_from(string.ascii_lowercase),
    )
    def test_combined_pattern_matches_either_alternative(
        self, name1: str, name2: str, ch1: str, ch2: str
    ) -> None:
        assume(name1 != name2)
        p1 = f"(?P<{name1}>[{ch1}])"
        p2 = f"(?P<{name2}>[{ch2}])"
        combined = re.compile(combine_patterns([p1, p2]))
        self.assertIsNotNone(combined.fullmatch(ch1))
        self.assertIsNotNone(combined.fullmatch(ch2))


if __name__ == "__main__":
    unittest.main()
