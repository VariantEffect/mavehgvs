import unittest
from mavehgvs.variant import Variant


class TestObjectCreation(unittest.TestCase):
    def test_protein_sub(self):
        variant_strings = ["p.Glu27Trp", "p.Ter345Lys", "p.Cys22="]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_dna_sub(self):
        variant_strings = ["g.48C>A", "c.122-6T>A", "c.*33G>C"]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_rna_sub(self):
        variant_strings = ["r.22g>u", "r.33+12a>c"]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_valid())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))

    def test_target_equivalent(self):
        variant_strings = [f"{prefix}.=" for prefix in "gmo" "cn" "r"]

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertTrue(v.is_target_identical())

        for s in variant_strings:
            with self.subTest(s=s):
                v = Variant(s)
                self.assertEqual(s, str(v))




if __name__ == '__main__':
    unittest.main()
