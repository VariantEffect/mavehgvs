import unittest
import re
from mavehgvs.patterns.dna import (
    dna_sub_cn,
    dna_sub_gmo,
    dna_del_cn,
    dna_del_gmo,
    dna_dup_cn,
    dna_dup_gmo,
    dna_ins_cn,
    dna_ins_gmo,
    dna_delins_cn,
    dna_delins_gmo,
    dna_variant_cn,
    dna_variant_gmo,
    dna_single_variant,
    dna_multi_variant,
)


class TestDnaSubCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_cn, flags=re.ASCII)

        cls.valid_strings = ["48C>A", "=", "122-6T>A", "*24G>C", "19+22A>G", "-27+3T>C"]

        cls.invalid_strings = ["22g>u", "48C>W", "22=", "122=/T>A"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaSubGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_gmo, flags=re.ASCII)

        cls.valid_strings = ["48C>A", "="]

        cls.invalid_strings = ["122-6T>A", "22g>u", "48C>W", "22=", "122=/T>A", "0C>T"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDelCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_del_cn, flags=re.ASCII)

        cls.valid_strings = [
            "44del",
            "1_95del",
            "78+5_78+10del",
            "-25+1_-25+3del",
            "*17del",
        ]

        cls.invalid_strings = [
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDelGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_del_gmo, flags=re.ASCII)

        cls.valid_strings = ["44del", "1_95del"]

        cls.invalid_strings = [
            "78+5_78+10del",
            "-25+1_-25+3del",
            "*17del",
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDupCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_dup_cn, flags=re.ASCII)

        cls.valid_strings = [
            "22_24dup",
            "77dup",
            "101+1_101+7dup",
            "-25+1_-25+3dup",
            "*17dup",
        ]

        cls.invalid_strings = [
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDupGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_dup_gmo, flags=re.ASCII)

        cls.valid_strings = ["22_24dup", "77dup"]

        cls.invalid_strings = [
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
            "101+1_101+7dup",
            "-25+1_-25+3dup",
            "*17dup",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaInsCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_ins_cn, flags=re.ASCII)

        cls.valid_strings = [
            "234_235insT",
            "84_85insCTG",
            "99+6_99+7insA",
            "124+100_124-100insTTG",
            "124+101_124-100insTTG",
        ]

        cls.invalid_strings = ["84_85ins100_125", "234_235ins(10)", "234_235ins(?)"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaInsGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_ins_gmo, flags=re.ASCII)

        cls.valid_strings = ["234_235insT", "84_85insCTG"]

        cls.invalid_strings = [
            "99+6_99+7insA",
            "84_85ins100_125",
            "234_235ins(10)",
            "234_235ins(?)",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDelinsCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_delins_cn, flags=re.ASCII)

        cls.valid_strings = [
            "22delinsAACG",
            "83_85delinsT",
            "43-6_595+12delinsCTT",
            "*788delinsA",
        ]

        cls.invalid_strings = ["84_85delinsAAN", "234delinsW"]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaDelinsGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_delins_gmo, flags=re.ASCII)

        cls.valid_strings = ["22delinsAACG", "83_85delinsT"]

        cls.invalid_strings = [
            "43-6_595+12delinsCTT",
            "*788delinsA",
            "84_85delinsAAN",
            "234delinsW",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaVariantCn(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_variant_cn, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "122-6T>A",
            "*24G>C",
            "19+22A>G",
            "-27+3T>C",
            "44del",
            "1_95del",
            "78+5_78+10del",
            "-25+1_-25+3del",
            "*17del",
            "22_24dup",
            "77dup",
            "101+1_101+7dup",
            "-25+1_-25+3dup",
            "*17dup",
            "234_235insT",
            "84_85insCTG",
            "99+6_99+7insA",
            "22delinsAACG",
            "83_85delinsT",
            "43-6_595+12delinsCTT",
            "*788delinsA",
        ]

        cls.invalid_strings = [
            "22g>u",
            "48C>W",
            "22=",
            "122=/T>A",
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
            "84_85ins100_125",
            "234_235ins(10)",
            "234_235ins(?)",
            "84_85delinsAAN",
            "234delinsW",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaVariantGmo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_variant_gmo, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "44del",
            "1_95del",
            "22_24dup",
            "77dup",
            "234_235insT",
            "84_85insCTG",
            "22delinsAACG",
            "83_85delinsT",
        ]

        cls.invalid_strings = [
            "43-6_595+12delinsCTT",
            "*788delinsA",
            "99+6_99+7insA",
            "101+1_101+7dup",
            "-25+1_-25+3dup",
            "*17dup",
            "78+5_78+10del",
            "-25+1_-25+3del",
            "*17del",
            "*24G>C",
            "19+22A>G",
            "122-6T>A",
            "-27+3T>C",
            "22g>u",
            "48C>W",
            "22=",
            "122=/T>A",
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
            "84_85ins100_125",
            "234_235ins(10)",
            "234_235ins(?)",
            "84_85delinsAAN",
            "234delinsW",
        ]

    def test_valid_strings(self):
        for s in self.valid_strings:
            with self.subTest(s=s):
                self.assertIsNotNone(
                    self.pattern.fullmatch(s), msg=f'failed to match "{s}"'
                )

    def test_invalid_strings(self):
        for s in self.invalid_strings:
            with self.subTest(s=s):
                self.assertIsNone(
                    self.pattern.fullmatch(s), msg=f'incorrectly matched "{s}"'
                )


class TestDnaSingleVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_single_variant, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "44del",
            "1_95del",
            "22_24dup",
            "77dup",
            "234_235insT",
            "84_85insCTG",
            "22delinsAACG",
            "83_85delinsT",
        ]

        cls.valid_strings_cn_only = [
            "43-6_595+12delinsCTT",
            "*788delinsA",
            "99+6_99+7insA",
            "101+1_101+7dup",
            "-25+1_-25+3dup",
            "*17dup",
            "78+5_78+10del",
            "-25+1_-25+3del",
            "*17del",
            "*24G>C",
            "19+22A>G",
            "122-6T>A",
            "-27+3T>C",
        ]

        cls.invalid_strings = [
            "22g>u",
            "48C>W",
            "22=",
            "122=/T>A",
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
            "84_85ins100_125",
            "234_235ins(10)",
            "234_235ins(?)",
            "84_85delinsAAN",
            "234delinsW",
        ]

    def test_valid_strings(self):
        for p in "cngmo":
            for s in self.valid_strings:
                with self.subTest(s=s, p=p):
                    v = f"{p}.{s}"
                    self.assertIsNotNone(
                        self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                    )
        for p in "cn":
            for s in self.valid_strings_cn_only:
                with self.subTest(s=s, p=p):
                    v = f"{p}.{s}"
                    self.assertIsNotNone(
                        self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                    )

    def test_invalid_strings(self):
        for p in "cngmo":
            for s in self.invalid_strings:
                with self.subTest(s=s, p=p):
                    v = f"{p}.{s}"
                    self.assertIsNone(
                        self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                    )
        for p in "gmo":
            for s in self.valid_strings_cn_only:
                with self.subTest(s=s, p=p):
                    v = f"{p}.{s}"
                    self.assertIsNone(
                        self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                    )


class TestDnaMultiVariant(unittest.TestCase):
    @unittest.expectedFailure
    def test_something(self):
        self.assertEqual(True, False)


if __name__ == "__main__":
    unittest.main()
