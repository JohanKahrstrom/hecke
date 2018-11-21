import unittest

from jk.hecke import permutation as p


class TestPermutation(unittest.TestCase):
    def testEquality(self):
        self.assertEqual(
            p.Permutation([3, 2, 1]),
            p.Permutation([3, 2, 1])
        )

        self.assertNotEqual(
            p.Permutation([3, 2, 1]),
            p.Permutation([3, 1, 2])
        )

        self.assertNotEqual(
            p.Permutation([3, 2, 1]),
            p.Permutation([2, 1])
        )

    def testComposition(self):
        result = (
            p.Permutation([2, 1, 3]) *
            p.Permutation([1, 3, 2])
        )

        self.assertEqual(
            result,
            p.Permutation([2, 3, 1])
        )

    def testLen(self):
        self.assertEqual(
            len(p.Permutation([3, 4, 2, 1])),
            4
        )

    def testHash(self):
        # Crappy test, but still
        self.assertNotEqual(
            hash(p.Permutation([1, 2, 3])),
            hash(p.Permutation([3, 2, 1]))
        )

    def testInverse(self):
        identity = p.Permutation([1, 2, 3])
        permutation = p.Permutation([1, 2, 3])
        self.assertEqual(permutation * permutation.inverse(), identity)

        permutation = p.Permutation([3, 2, 1])
        self.assertEqual(permutation * permutation.inverse(), identity)

        permutation = p.Permutation([2, 3, 1])
        self.assertEqual(permutation * permutation.inverse(), identity)

        identity = p.Permutation([1, 2, 3, 4, 5, 6, 7, 8])
        permutation = p.Permutation([5, 6, 3, 1, 2, 8, 4, 7])
        self.assertEqual(permutation * permutation.inverse(), identity)
