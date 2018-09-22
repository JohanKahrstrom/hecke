import unittest
import hecke.hecke as h
import hecke.coxeter as c
import hecke.laurent as l
import hecke.permutation as p

class TestHecke(unittest.TestCase):
    def test_add(self):
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

    def test_mul_identity(self):
        group = c.generate_a2()

        e = h.Hecke(group, {
            'e': l.Laurent({0: 1})
        })

        self.assertEqual(e * e, e)

        h1 = h.Hecke(group, {
            's': l.Laurent({10: 1, 20: 2}),
            'rs': l.Laurent({-1: 2, 0: 1})
        })

        self.assertEqual(h1 * e, h1)

        h2 = h.Hecke(group, {
            'e': l.Laurent({-100: 20, 0: 10, 100: 20}),
            'rsr': l.Laurent({-25: -12, -2: -1, 2: 2})
        })

        self.assertEqual(h2 * e, h2)

    def test_mul_generator(self):
        group = c.generate_a2()

        e = h.Hecke(group, {
            'e': l.Laurent({0: 1})
        })
        r = h.Hecke(group, {
            'r': l.Laurent({0: 1})
        })
        s = h.Hecke(group, {
            's': l.Laurent({0: 1})
        })
        rs = h.Hecke(group, {
            'rs': l.Laurent({0: 1})
        })
        sr = h.Hecke(group, {
            'sr': l.Laurent({0: 1})
        })
        rsr = h.Hecke(group, {
            'rsr': l.Laurent({0: 1})
        })


#        self.assertEqual(e * s, s)

#        self.assertEqual(e * r, r)

#        self.assertEqual(r * s, rs)

#        self.assertEqual(s * r, sr)

#        self.assertEqual(r * s * r, rsr)
#        self.assertEqual(rs * r, rsr)
#        self.assertEqual(r * sr, rsr)
#        self.assertEqual(sr * s, rsr)