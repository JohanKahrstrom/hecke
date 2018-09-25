import unittest

import hecke.permutation as p
import hecke.coxeter as c


class TestCoxeter(unittest.TestCase):
    def testGenerateFailsOnBadInput(self):
        self.assertRaises(Exception,
                          c.CoxeterGroup.generate, {
                              'a': p.Permutation([1]),
                              'b': p.Permutation([2, 1])
                          })

    def testGenerateReturnsIdentity(self):
        group = c.generate_a1()

        self.assertTrue('e' in group.elements)
        self.assertEquals(group.name_lookup[p.Permutation([1, 2])], 'e')

    def testGeneratesAllElementsA1(self):
        group = c.generate_a1()

        expected_lookup = {
            p.Permutation([1, 2]): 'e',
            p.Permutation([2, 1]): 's'
        }

        self.assertEqual(group.name_lookup, expected_lookup)

    def testMultiplicationInA1(self):
        group = c.generate_a1()
        s = group['s']
        ss = s * s
        self.assertEqual(ss, group['e'])
        self.assertNotEqual(s, ss)

    def testLength(self):
        group = c.generate_a2()

        self.assertEqual(group['e'].length(), 0)
        self.assertEqual(group['r'].length(), 1)
        self.assertEqual(group['s'].length(), 1)
        self.assertEqual(group['rs'].length(), 2)
        self.assertEqual(group['sr'].length(), 2)
        self.assertEqual(group['rsr'].length(), 3)

    def testGeneratesAllElementsA2(self):
        group = c.generate_a2()

        expected_lookup = {
            p.Permutation([1, 2, 3]): 'e',
            p.Permutation([2, 1, 3]): 'r',
            p.Permutation([1, 3, 2]): 's',
            p.Permutation([2, 3, 1]): 'rs',
            p.Permutation([3, 1, 2]): 'sr',
            p.Permutation([3, 2, 1]): 'rsr'
        }

        self.assertEqual(group.name_lookup, expected_lookup)

        e = group['e']
        r = group['r']
        s = group['s']
        rs = group['rs']
        sr = group['sr']
        rsr = group['rsr']

        self.assertEqual(r*r, e)
        self.assertEqual(s*s, e)
        self.assertEqual(r*s, rs)
        self.assertEqual(s*r, sr)
        self.assertEqual(sr*r, s)
        self.assertEqual(rs*s, r)
        self.assertEqual(s*sr, r)
        self.assertEqual(r*rs, s)
        self.assertEqual(r*s*r, rsr)
        self.assertEqual(s*r*s, rsr)
        self.assertEqual(rs*r, rsr)
        self.assertEqual(r*sr, rsr)
        self.assertEqual(s*rs, rsr)
        self.assertEqual(sr*s, rsr)

    def testInverse(self):
        group = c.generate_a2()

        for element in group.elements.values():
            self.assertEqual(element * element.inverse(), group.identity)
            self.assertEqual(element.inverse() * element, group.identity)

    def test_longest(self):
        a1 = c.generate_a1()
        self.assertEqual(a1.longest.name, "s")

        a2 = c.generate_a2()
        self.assertEqual(a2.longest.name, "rsr")

        a3 = c.generate_a3()
        self.assertEqual(a3.longest.name, "rsrtsr")