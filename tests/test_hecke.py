import unittest
import hecke.hecke as h
import hecke.coxeter as c
import hecke.laurent as l


class TestHecke(unittest.TestCase):
    def test_add(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        h1 = h.HeckeElement(hecke, {
            's': l.Laurent({-1: 1, 1: 1})
        })
        h2 = h.HeckeElement(hecke, {
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({0: 1, 1: 4})
        })

        expected_sum = h.HeckeElement(hecke, {
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({-1: 1, 0: 1, 1: 5})
        })

        self.assertEqual(h1 + h2, expected_sum)

    def test_mul_identity(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        e = h.HeckeElement(hecke, {
            'e': l.Laurent({0: 1})
        })

        self.assertEqual(e * e, e)

        h1 = h.HeckeElement(hecke, {
            's': l.Laurent({10: 1, 20: 2}),
            'rs': l.Laurent({-1: 2, 0: 1})
        })

        self.assertEqual(h1 * e, h1)

        h2 = h.HeckeElement(hecke, {
            'e': l.Laurent({-100: 20, 0: 10, 100: 20}),
            'rsr': l.Laurent({-25: -12, -2: -1, 2: 2})
        })

        self.assertEqual(h2 * e, h2)
        self.assertEqual(e * h2, h2)

    def test_mul_generator(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        e = h.HeckeElement(hecke, {
            'e': l.Laurent({0: 1})
        })
        r = h.HeckeElement(hecke, {
            'r': l.Laurent({0: 1})
        })
        s = h.HeckeElement(hecke, {
            's': l.Laurent({0: 1})
        })
        rs = h.HeckeElement(hecke, {
            'rs': l.Laurent({0: 1})
        })
        sr = h.HeckeElement(hecke, {
            'sr': l.Laurent({0: 1})
        })
        rsr = h.HeckeElement(hecke, {
            'rsr': l.Laurent({0: 1})
        })

        self.assertEqual(e * s, s)
        self.assertEqual(e * r, r)
        self.assertEqual(r * s, rs)
        self.assertEqual(s * r, sr)
        self.assertEqual(r * s * r, rsr)
        self.assertEqual(rs * r, rsr)
        self.assertEqual(r * sr, rsr)
        self.assertEqual(sr * s, rsr)
        self.assertEqual(s * s,
                         h.HeckeElement(hecke, {
                             's': l.Laurent({-1: 1, 1: -1}),
                             'e': l.Laurent({0: 1})
                         }))
        self.assertEqual(sr * r,
                         h.HeckeElement(hecke, {
                             'sr': l.Laurent({-1: 1, 1: -1}),
                             's': l.Laurent({0: 1})
                         }))
        self.assertEqual(h.HeckeElement(hecke, {
                             'r': l.Laurent({0: 2})
                         }) * h.HeckeElement(group, {
                             'sr': l.Laurent({0: 2})
                         }),
                         h.HeckeElement(hecke, {
                             'rsr': l.Laurent({0: 4})
                         }))

    def test_w0_central(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)
        w0 = h.HeckeElement(hecke, {'rsr': l.Laurent({0: 1})})
        w02 = w0 * w0
        h1 = h.HeckeElement(hecke, {'e': l.Laurent({1: 2}),
                                    'sr': l.Laurent({-10: 20, 2: 1})})
        h2 = h.HeckeElement(hecke, {'e': l.Laurent({0: 1})})

        self.assertEqual(w02 * h1, h1 * w02)
        self.assertEqual(w02 * h2, h2 * w02)

    def test_get_standard_base(self):
        group = c.generate_a2()

        for name, element in group.elements.items():
            pass