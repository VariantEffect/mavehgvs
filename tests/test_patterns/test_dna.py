import unittest
import re
from hypothesis import given
from mavehgvs.patterns.dna import (
    dna_equal_c,
    dna_equal_n,
    dna_equal_gmo,
    dna_sub_c,
    dna_sub_n,
    dna_sub_gmo,
    dna_del_c,
    dna_del_n,
    dna_del_gmo,
    dna_dup_c,
    dna_dup_n,
    dna_dup_gmo,
    dna_ins_c,
    dna_ins_n,
    dna_ins_gmo,
    dna_delins_c,
    dna_delins_n,
    dna_delins_gmo,
    dna_variant_c,
    dna_variant_n,
    dna_variant_gmo,
    dna_single_variant,
    dna_multi_variant,
)
from . import (
    build_multi_variants,
    plain_positions,
    intron_offset_positions,
    intron_only_positions,
    utr_intron_positions,
    utr_only_positions,
    nucleotide_variant_strings,
)


class TestDnaEqualC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_equal_c, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "18=",
            "10_14=",
            "122-6=",
            "*24=",
            "19+22=",
            "19+22_88=",
            "-27+3=",
        ]

        cls.invalid_strings = ["=22", "(=)", "18(=)"]

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


class TestDnaEqualN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_equal_n, flags=re.ASCII)

        cls.valid_strings = [
            "=",
            "18=",
            "10_14=",
            "122-6=",
            "19+22=",
            "19+22_88=",
        ]

        cls.invalid_strings = [
            "=22",
            "(=)",
            "18(=)",
            "-27+3=",
            "*24=",
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


class TestDnaEqualGMO(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_equal_gmo, flags=re.ASCII)

        cls.valid_strings = ["=", "18=", "10_14="]

        cls.invalid_strings = [
            "=22",
            "(=)",
            "18(=)",
            "122-6=",
            "*24=",
            "19+22=",
            "19+22_88=",
            "-27+3=",
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


class TestDnaSubC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_c, flags=re.ASCII)

        cls.valid_strings = ["48C>A", "122-6T>A", "*24G>C", "19+22A>G", "-27+3T>C"]

        cls.invalid_strings = ["22g>u", "48C>W", "122=/T>A"]

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


class TestDnaSubN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_sub_n, flags=re.ASCII)

        cls.valid_strings = ["48C>A", "122-6T>A", "19+22A>G"]

        cls.invalid_strings = ["22g>u", "48C>W", "122=/T>A", "*24G>C", "-27+3T>C"]

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

        cls.valid_strings = ["48C>A"]

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


class TestDnaDelC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_del_c, flags=re.ASCII)

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


class TestDnaDelN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_del_n, flags=re.ASCII)

        cls.valid_strings = ["44del", "1_95del", "78+5_78+10del"]

        cls.invalid_strings = [
            "(78+1_79-1)_(124+1_125-1)del",
            "(?_85)_(124_?)del",
            "122=/del",
            "-25+1_-25+3del",
            "*17del",
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


class TestDnaDupC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_dup_c, flags=re.ASCII)

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


class TestDnaDupN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_dup_n, flags=re.ASCII)

        cls.valid_strings = ["22_24dup", "77dup", "101+1_101+7dup"]

        cls.invalid_strings = [
            "(78+1_79-1)_(124+1_125-1)dup",
            "(?_85)_(124_?)dup",
            "122_125=//dup",
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


class TestDnaInsC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_ins_c, flags=re.ASCII)

        cls.valid_strings = [
            "234_235insT",
            "84_85insCTG",
            "*84_*85insCTG",
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


class TestDnaInsN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_ins_n, flags=re.ASCII)

        cls.valid_strings = [
            "234_235insT",
            "84_85insCTG",
            "99+6_99+7insA",
            "124+100_124-100insTTG",
            "124+101_124-100insTTG",
        ]

        cls.invalid_strings = [
            "84_85ins100_125",
            "234_235ins(10)",
            "234_235ins(?)",
            "*84_*85insCTG",
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


class TestDnaDelinsC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_delins_c, flags=re.ASCII)

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


class TestDnaDelinsN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_delins_n, flags=re.ASCII)

        cls.valid_strings = ["22delinsAACG", "83_85delinsT", "43-6_595+12delinsCTT"]

        cls.invalid_strings = ["84_85delinsAAN", "234delinsW" "*788delinsA"]

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


class TestDnaVariantC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_variant_c, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "22=",
            "4_6=",
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


class TestDnaVariantN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_variant_n, flags=re.ASCII)

        cls.valid_strings = [
            "48C>A",
            "=",
            "22=",
            "1_3=",
            "122-6T>A",
            "19+22A>G",
            "44del",
            "1_95del",
            "78+5_78+10del",
            "22_24dup",
            "77dup",
            "101+1_101+7dup",
            "234_235insT",
            "84_85insCTG",
            "99+6_99+7insA",
            "22delinsAACG",
            "83_85delinsT",
            "43-6_595+12delinsCTT",
        ]

        cls.invalid_strings = [
            "22g>u",
            "48C>W",
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
            "*24G>C",
            "-27+3T>C",
            "-25+1_-25+3del",
            "*17del",
            "-25+1_-25+3dup",
            "*17dup",
            "*788delinsA",
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
            "22=",
            "1_3=",
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
            "22=",
            "4_6=",
            "44del",
            "1_95del",
            "22_24dup",
            "77dup",
            "234_235insT",
            "84_85insCTG",
            "22delinsAACG",
            "83_85delinsT",
        ]

        cls.valid_strings_c_only = [
            "*788delinsA",
            "-25+1_-25+3dup",
            "*17dup",
            "-25+1_-25+3del",
            "*17del",
            "*24G>C",
            "-27+3T>C",
        ]

        cls.valid_strings_cn_only = [
            "43-6_595+12delinsCTT",
            "99+6_99+7insA",
            "101+1_101+7dup",
            "78+5_78+10del",
            "19+22A>G",
            "122-6T>A",
        ]

        cls.invalid_strings = [
            "22g>u",
            "48C>W",
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
        for p in "c":
            for s in self.valid_strings_c_only:
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
        for p in "ngmo":
            for s in self.valid_strings_c_only:
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
    @classmethod
    def setUpClass(cls):
        cls.pattern = re.compile(dna_multi_variant, flags=re.ASCII)

        single_valid_strings = [
            "48C>A",
            "=",
            "22=",
            "4_6=",
            "44del",
            "1_95del",
            "22_24dup",
            "77dup",
            "234_235insT",
            "84_85insCTG",
            "22delinsAACG",
            "83_85delinsT",
        ]

        single_valid_strings_c_only = [
            "*788delinsA",
            "-25+1_-25+3dup",
            "*17dup",
            "-25+1_-25+3del",
            "*17del",
            "*24G>C",
            "-27+3T>C",
        ]

        single_valid_strings_cn_only = [
            "43-6_595+12delinsCTT",
            "99+6_99+7insA",
            "101+1_101+7dup",
            "78+5_78+10del",
            "19+22A>G",
            "122-6T>A",
        ]

        single_invalid_strings = [
            "22g>u",
            "48C>W",
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

        cls.valid_strings, cls.invalid_strings = build_multi_variants(
            single_valid_strings, single_invalid_strings
        )
        cls.valid_strings_c_only, cls.invalid_strings_ngmo = build_multi_variants(
            single_valid_strings_c_only, single_valid_strings_c_only
        )
        cls.valid_strings_cn_only, cls.invalid_strings_gmo = build_multi_variants(
            single_valid_strings_cn_only, single_valid_strings_cn_only
        )

    def test_valid_strings(self):
        for p in "cngmo":
            for s in self.valid_strings:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNotNone(
                        self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                    )
        for p in "cn":
            for s in self.valid_strings_cn_only:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNotNone(
                        self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                    )
        for p in "c":
            for s in self.valid_strings_c_only:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNotNone(
                        self.pattern.fullmatch(v), msg=f'failed to match "{v}"'
                    )

    def test_invalid_strings(self):
        for p in "cngmo":
            for s in self.invalid_strings:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNone(
                        self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                    )
        for p in "ngmo":
            for s in self.invalid_strings_ngmo:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNone(
                        self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                    )
        for p in "gmo":
            for s in self.invalid_strings_gmo:
                with self.subTest(s=s, p=p):
                    v = f"{p}.[{s}]"
                    self.assertIsNone(
                        self.pattern.fullmatch(v), msg=f'incorrectly matched "{v}"'
                    )


class TestDnaHypothesis(unittest.TestCase):
    """Property-based tests that generalize the fixed examples above: any
    generated variant string for a given position flavor should match the
    corresponding combined variant pattern and be accepted by dna_single_variant
    under the correct prefix(es); positions using features outside a flavor's
    grammar (UTR / intron) should be rejected by the narrower flavors.
    """

    @classmethod
    def setUpClass(cls):
        cls.variant_c = re.compile(dna_variant_c, flags=re.ASCII)
        cls.variant_n = re.compile(dna_variant_n, flags=re.ASCII)
        cls.variant_gmo = re.compile(dna_variant_gmo, flags=re.ASCII)
        cls.single = re.compile(dna_single_variant, flags=re.ASCII)
        cls.multi = re.compile(dna_multi_variant, flags=re.ASCII)

    @given(s=nucleotide_variant_strings(utr_intron_positions(), "ACGT"))
    def test_generated_c_variants(self, s: str) -> None:
        self.assertIsNotNone(self.variant_c.fullmatch(s), msg=f'failed to match "{s}"')
        self.assertIsNotNone(
            self.single.fullmatch(f"c.{s}"), msg=f'failed to match "c.{s}"'
        )

    @given(s=nucleotide_variant_strings(intron_offset_positions(), "ACGT"))
    def test_generated_n_variants(self, s: str) -> None:
        self.assertIsNotNone(self.variant_n.fullmatch(s), msg=f'failed to match "{s}"')
        self.assertIsNotNone(
            self.single.fullmatch(f"n.{s}"), msg=f'failed to match "n.{s}"'
        )

    @given(s=nucleotide_variant_strings(plain_positions(), "ACGT"))
    def test_generated_gmo_variants(self, s: str) -> None:
        self.assertIsNotNone(
            self.variant_gmo.fullmatch(s), msg=f'failed to match "{s}"'
        )
        for p in "gmo":
            self.assertIsNotNone(
                self.single.fullmatch(f"{p}.{s}"), msg=f'failed to match "{p}.{s}"'
            )

    @given(
        s=nucleotide_variant_strings(
            utr_only_positions(), "ACGT", allow_bare_equal=False
        )
    )
    def test_utr_positions_rejected_outside_c(self, s: str) -> None:
        self.assertIsNone(self.variant_n.fullmatch(s), msg=f'incorrectly matched "{s}"')
        self.assertIsNone(
            self.variant_gmo.fullmatch(s), msg=f'incorrectly matched "{s}"'
        )

    @given(
        s=nucleotide_variant_strings(
            intron_only_positions(), "ACGT", allow_bare_equal=False
        )
    )
    def test_intron_positions_rejected_by_gmo(self, s: str) -> None:
        self.assertIsNone(
            self.variant_gmo.fullmatch(s), msg=f'incorrectly matched "{s}"'
        )

    @given(
        s1=nucleotide_variant_strings(utr_intron_positions(), "ACGT"),
        s2=nucleotide_variant_strings(utr_intron_positions(), "ACGT"),
    )
    def test_generated_multi_variant(self, s1: str, s2: str) -> None:
        v = f"c.[{s1};{s2}]"
        self.assertIsNotNone(self.multi.fullmatch(v), msg=f'failed to match "{v}"')


if __name__ == "__main__":
    unittest.main()
