import unittest
import hecke.hecke as h
import hecke.coxeter as c
import hecke.laurent as l
import hecke.permutation as p

class TestHecke(unittest.TestCase):
    def testAdd(self):
        group = c.generate_a2()

        h1 = h.Hecke(group, {
            's': l.Laurent({-1: 1, 1: 1})
        })
        h2 = h.Hecke(group, {
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({0: 1, 1: 4})
        })

        expected_sum = h.Hecke(group, {
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({-1: 1, 0: 1, 1: 5})
        })

        self.assertEqual(h1 + h2, expected_sum)