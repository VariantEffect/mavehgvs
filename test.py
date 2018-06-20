import unittest

from hgvsp.tests import test_hgvs, test_dna_hgvs, \
    test_protein_hgvs, test_rna_hgvs


if __name__ == "__main__":
    loader = unittest.TestLoader()
    tests = loader.discover(start_dir='./', pattern="test_*.py")
    unittest.TextTestRunner().run(tests)
