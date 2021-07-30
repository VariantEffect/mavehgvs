import unittest

from mavehgvs.patterns.util import combine_patterns, remove_named_groups


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


if __name__ == "__main__":
    unittest.main()
