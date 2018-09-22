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
        group = c.CoxeterGroup.generate({'s': p.Permutation([2, 1])})

        self.assertTrue('e' in group.elements)
        self.assertEquals(group.name_lookup[p.Permutation([1, 2])], 'e')

    def testGeneratesAllElementsA1(self):
        group = c.CoxeterGroup.generate({'s': p.Permutation([2, 1])})

        expected_lookup = {
            p.Permutation([1, 2]): 'e',
            p.Permutation([2, 1]): 's'
        }

        self.assertEqual(group.name_lookup, expected_lookup)

    def testMultiplicationInA1(self):
        group = c.CoxeterGroup.generate({'s': p.Permutation([2, 1])})
        s = group.get('s')
        ss = s * s
        self.assertEqual(ss, group.get('e'))
        self.assertNotEqual(s, ss)

    def testLength(self):
        group = c.generate_a2()

        self.assertEqual(group.get('e').length(), 0)
        self.assertEqual(group.get('r').length(), 1)
        self.assertEqual(group.get('s').length(), 1)
        self.assertEqual(group.get('rs').length(), 2)
        self.assertEqual(group.get('sr').length(), 2)
        self.assertEqual(group.get('rsr').length(), 3)

    def testGeneratesAllElementsA2(self):
        group = c.CoxeterGroup.generate({
            'r': p.Permutation([2, 1, 3]),
            's': p.Permutation([1, 3, 2])
        })

        expected_lookup = {
            p.Permutation([1, 2, 3]): 'e',
            p.Permutation([2, 1, 3]): 'r',
            p.Permutation([1, 3, 2]): 's',
            p.Permutation([2, 3, 1]): 'rs',
            p.Permutation([3, 1, 2]): 'sr',
            p.Permutation([3, 2, 1]): 'rsr'
        }

        self.assertEqual(group.name_lookup, expected_lookup)

        e = group.get('e')
        r = group.get('r')
        s = group.get('s')
        rs = group.get('rs')
        sr = group.get('sr')
        rsr = group.get('rsr')

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
