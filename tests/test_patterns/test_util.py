import unittest


class TestCombinePatterns(unittest.TestCase):
    @unittest.expectedFailure
    def test_something(self):
        self.assertEqual(True, False)


class TestRemoveNamedGroups(unittest.TestCase):
    @unittest.expectedFailure
    def test_something(self):
        self.assertEqual(True, False)


if __name__ == "__main__":
    unittest.main()
