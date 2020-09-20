import unittest
import re
from mavehgvs.patterns.dna import dna_sub_gmo, dna_sub_cn


class TestDnaSubCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_cn, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "122-6T>A",
            "*24G>C",
            "19+22A>G",
        ]

        cls.invalid_strings = [
            "22g>u",
            "48C>W",
            "22=",
            "122=/T>A",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(self.pattern.fullmatch(s))

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(self.pattern.fullmatch(s))


class TestDnaSubGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_gmo, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
        ]

        cls.invalid_strings = [
            "122-6T>A",
            "22g>u",
            "48C>W",
            "22=",
            "122=/T>A",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(self.pattern.fullmatch(s))

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(self.pattern.fullmatch(s))


if __name__ == "__main__":
    unittest.main()
