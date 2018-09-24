import unittest

import hecke.laurent as l


class TestLaurent(unittest.TestCase):
    def testEquals(self):
        self.assertEqual(
            l.Laurent({0: 1, 1: 2}),
            l.Laurent({0: 1, 1: 2})
        )
        self.assertNotEqual(
            l.Laurent({0: 1, 1: 2}),
            l.Laurent({0: 1, 1: 3})
        )

    def testAdd(self):
        result = (
            l.Laurent({0: 10, 1: 2}) +
            l.Laurent({-1: 1, 0: 2, 1: 3, 2: 3})
        )

        self.assertEqual(
            result,
            l.Laurent({-1: 1, 0: 12, 1: 5, 2: 3})
        )

    def testNeg(self):
        result = -l.Laurent({-1: 1, 0: -10, 2: 15})

        self.assertEqual(
            result,
            l.Laurent({-1: -1, 0: 10, 2: -15})
        )

    def testSub(self):
        result = (
            l.Laurent({-1: 1, 0: 2, 1: 3}) -
            l.Laurent({-2: 10, -1: -1, 1: 2})
        )

        self.assertEqual(
            result,
            l.Laurent({-2: -10, -1: 2, 0: 2, 1: 1})
        )

    def testMult(self):
        result = (
            l.Laurent({-1: 1, 0: 2, 1: 3}) *
            l.Laurent({0: 2, 1: 4, 2: -1})
        )

        self.assertEqual(
            result,
            l.Laurent({
                -1: 1*2,
                0: 2*2 + 1*4,
                1: 1*(-1) + 2*4 + 3*2,
                2: 2*(-1) + 3*4,
                3: 3*(-1)
            })
        )

    def testMultInt(self):
        result = l.Laurent({-4: 2, -1: 0, 2: 3}) * 3

        self.assertEqual(
            result,
            l.Laurent({
                -4: 6,
                -1: 0,
                2: 9
            })
        )

    def testInvolution(self):
        p = l.Laurent({-4: 2, -1: 0, 2: 3})
        result = p.involute()
        self.assertEqual(
            result,
            l.Laurent({-2: 3, 1: 0, 4: 2})
        )

        self.assertEqual(
            p.involute().involute(),
            p
        )