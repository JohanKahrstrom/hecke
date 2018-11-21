import unittest
from jk.hecke import hecke as h, coxeter as c
import jk.hecke.laurent as l
import numpy as np


class TestHecke(unittest.TestCase):
    def test_add(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        h1 = hecke.element({
            's': l.Laurent({-1: 1, 1: 1})
        })
        h2 = hecke.element({
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({0: 1, 1: 4})
        })

        expected_sum = hecke.element({
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({-1: 1, 0: 1, 1: 5})
        })

        self.assertEqual(h1 + h2, expected_sum)

    def test_sub(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        h1 = hecke.element({
            's': l.Laurent({-1: 1, 1: 1})
        })

        h2 = hecke.element({
            'e': l.Laurent({-3: 5}),
            's': l.Laurent({0: 1, 1: 4})
        })

        h3 = hecke.element(({
            'e': l.Laurent({-3: -5}),
            's': l.Laurent({-1: 1, 0: -1, 1: -3})
        }))

        self.assertEqual(h1 - h1, hecke.zero)
        self.assertEqual(h2 - h2, hecke.zero)

        self.assertEqual(h1 - h2, h3)

    def test_mul_identity(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        e = hecke.element({
            'e': l.Laurent({0: 1})
        })

        self.assertEqual(e * e, e)

        h1 = hecke.element({
            's': l.Laurent({10: 1, 20: 2}),
            'rs': l.Laurent({-1: 2, 0: 1})
        })

        self.assertEqual(h1 * e, h1)

        h2 = hecke.element({
            'e': l.Laurent({-100: 20, 0: 10, 100: 20}),
            'rsr': l.Laurent({-25: -12, -2: -1, 2: 2})
        })

        self.assertEqual(h2 * e, h2)
        self.assertEqual(e * h2, h2)

    def test_mul_generator(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        e = hecke.element({
            'e': l.Laurent({0: 1})
        })
        r = hecke.element({
            'r': l.Laurent({0: 1})
        })
        s = hecke.element({
            's': l.Laurent({0: 1})
        })
        rs = hecke.element({
            'rs': l.Laurent({0: 1})
        })
        sr = hecke.element({
            'sr': l.Laurent({0: 1})
        })
        rsr = hecke.element({
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
                         hecke.element({
                             's': l.Laurent({-1: 1, 1: -1}),
                             'e': l.Laurent({0: 1})
                         }))
        self.assertEqual(sr * r,
                         hecke.element({
                             'sr': l.Laurent({-1: 1, 1: -1}),
                             's': l.Laurent({0: 1})
                         }))
        self.assertEqual(hecke.element({
                             'r': l.Laurent({0: 2})
                         }) * hecke.element({
                             'sr': l.Laurent({0: 2})
                         }),
                         hecke.element({
                             'rsr': l.Laurent({0: 4})
                         }))

    def test_mult_int(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        element = hecke.element({
            'e': l.Laurent({-100: 20, 0: 10, 100: 20}),
            'rsr': l.Laurent({-25: -12, -2: -1, 2: 2})
        })
        result = element * 3
        expected_result = hecke.element({
            'e': l.Laurent({-100: 60, 0: 30, 100: 60}),
            'rsr': l.Laurent({-25: -36, -2: -3, 2: 6})
        })

        self.assertEqual(result, expected_result)

    def test_mult_laurent(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        l1 = l.Laurent({0: 1, 1: -1})
        l2 = l.Laurent({-2: 10, 3: 5})

        element = hecke.element({
            'e': l1
        })
        result = element * l2
        expected_result = hecke.element({
            'e': l1 * l2
        })
        self.assertEqual(result,
                         expected_result)

    def test_w0_central(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)
        w0 = hecke.element({'rsr': l.Laurent({0: 1})})
        w02 = w0 * w0
        h1 = hecke.element({'e': l.Laurent({1: 2}),
                            'sr': l.Laurent({-10: 20, 2: 1})})
        h2 = hecke.element({'e': l.Laurent({0: 1})})

        self.assertEqual(w02 * h1, h1 * w02)
        self.assertEqual(w02 * h2, h2 * w02)

    def test_get_inverse(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        self.assertEqual(
            hecke.get_standard_basis_element('e') * hecke.get_generator_inverse_element('e'),
            hecke.one
        )

        self.assertEqual(
            hecke.get_standard_basis_element('r') * hecke.get_generator_inverse_element('r'),
            hecke.one
        )

        self.assertEqual(
            hecke.get_standard_basis_element('s') * hecke.get_generator_inverse_element('s'),
            hecke.one
        )

        self.assertEqual(
            hecke.get_standard_basis_element('rs') * hecke.get_standard_inverse_element('rs'),
            hecke.one
        )

        self.assertEqual(
            hecke.get_standard_basis_element('sr') * hecke.get_standard_inverse_element('sr'),
            hecke.one
        )

        self.assertEqual(
            hecke.get_standard_basis_element('rsr') * hecke.get_standard_inverse_element('rsr'),
            hecke.one
        )

    def test_dual(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        self.assertEqual(
            hecke.get_standard_basis_element('e').dual().dual(),
            hecke.get_standard_basis_element('e')
        )
        self.assertEqual(
            hecke.get_standard_basis_element('r').dual().dual(),
            hecke.get_standard_basis_element('r')
        )
        self.assertEqual(
            hecke.get_standard_basis_element('s').dual().dual(),
            hecke.get_standard_basis_element('s')
        )
        self.assertEqual(
            hecke.get_standard_basis_element('rs').dual().dual(),
            hecke.get_standard_basis_element('rs')
        )
        self.assertEqual(
            hecke.get_standard_basis_element('sr').dual().dual(),
            hecke.get_standard_basis_element('sr')
        )
        self.assertEqual(
            hecke.get_standard_basis_element('rsr').dual().dual(),
            hecke.get_standard_basis_element('rsr')
        )

        h1 = hecke.element({
            'e': l.Laurent({-1: 10, 2: 3}),
            's': l.Laurent({0: 1, 1: 3}),
            'rsr': l.Laurent({20: -12})
        })

        self.assertEqual(h1.dual().dual(), h1)

        hs = hecke.element({
            'e': l.Laurent({1: 1}),
            's': l.Laurent({0: 1})
        })

        self.assertEqual(hs.dual(), hs)
        print(hs * hs)
        self.assertEqual(
            hs * hs,
            hs * l.Laurent({-1: 1, 1: 1})
        )

    def assert_positive(self, element, head):
        self.assertEqual(element[head], l.one)
        for key, value in element.elements.items():
            if key != head:
                self.assertTrue(value.all_positive_degree())

    def test_generate_kl_basis(self):

        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        kl_basis = hecke.generate_kl_basis()

        for name in group.elements.keys():
            self.assertEqual(kl_basis[name], kl_basis[name].dual())
            self.assert_positive(kl_basis[name], name)

    def test_tau(self):
        group = c.generate_a3()
        hecke = h.HeckeAlgebra(group)

        for x in group.elements.values():
            hx = hecke.get_standard_basis_element(x.name)
            for y in group.elements.values():
                hy = hecke.get_standard_basis_element(y.name)
                if x == y.inverse():
                    self.assertEqual(
                        (hx * hy).tau(),
                        l.one
                    )
                else:
                    self.assertEqual(
                        (hx * hy).tau(),
                        l.zero
                    )

    def test_generate_dual_kl_basis(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        kl_basis = hecke.generate_kl_basis()
        dual_kl_basis = hecke.generate_dual_kl_basis()

        for x in group.elements.values():
            kl_x = kl_basis[x.name]
            for y in group.elements.values():
                dkl_y = dual_kl_basis[y.name]
                if x == y.inverse():
                    self.assertEqual(
                        (kl_x * dkl_y).tau(),
                        l.one
                    )
                    self.assertEqual(
                        (dkl_y * kl_x).tau(),
                        l.one
                    )
                else:
                    self.assertEqual(
                        (kl_x * dkl_y).tau(),
                        l.zero
                    )
                    self.assertEqual(
                        (dkl_y * kl_x).tau(),
                        l.zero
                    )

        for name in group.elements.keys():
            self.assertEqual(kl_basis[name], kl_basis[name].dual())
            self.assert_positive(kl_basis[name], name)

    def test_i(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        for element in group.all_elements():
            self.assertEqual(
                hecke.get_standard_basis_element(element.name).i(),
                hecke.get_standard_basis_element(element.inverse().name)
            )

    def test_kl_basis(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)
        kl_basis = hecke.generate_kl_basis()

        for x in group.elements.keys():
            hx = hecke.get_standard_basis_element(x)
            hxkl = hx.in_kl_basis()
            ress = hecke.zero
            for element, coef in hxkl.items():
                ress += kl_basis[element] * coef
            self.assertEqual(hx, ress)

    def test_dual_kl_basis(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)
        dual_kl_basis = hecke.generate_dual_kl_basis()

        for x in group.elements.keys():
            hx = hecke.get_standard_basis_element(x)
            hxkl = hx.in_dual_kl_basis()
            ress = hecke.zero
            for element, coef in hxkl.items():
                ress += dual_kl_basis[element] * coef
            self.assertEqual(hx, ress)

    def test_filtration_string(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)
        kl_basis = hecke.generate_kl_basis()
        dual_kl_basis = hecke.generate_dual_kl_basis()

        self.assertEqual(
            dual_kl_basis['e'].dual_kl_filtration(),
            'e\n'
        )

        self.assertEqual(
            hecke.get_standard_basis_element('e').dual_kl_filtration(),
            (
                "  e  \n"
                " r s \n"
                "rs sr\n"
                " rsr \n"
            )
        )

    def test_basis_matrix(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        standard_basis = [
            (x, hecke.get_standard_basis_element(x))
            for x in sorted(group.elements.keys(), key=lambda x: (len(x), x))
        ]
        ret = hecke.basis_matrix(standard_basis,
                                 hecke.generate_dual_kl_basis())

        expected = [
            [1, 1, 1, 1, 1, 1],
            [0, 1, 0, 1, 1, 1],
            [0, 0, 1, 1, 1, 1],
            [0, 0, 0, 1, 0, 1],
            [0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 1]
        ]
        self.assertEquals(ret, expected)

        group = c.generate_a3()
        hecke = h.HeckeAlgebra(group)

        kl_basis = hecke.generate_kl_basis()
        dual_kl_basis = hecke.generate_dual_kl_basis()

        d = dict()
        for element in ['e', 's', 'sr', 'st']:
            d[element] = (dual_kl_basis['s'] * kl_basis[element]).in_dual_kl_basis()

        ret = hecke.matrix(d)

        expected = [
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]

        self.assertEqual(ret, expected)

    def test_orders(self):
        group = c.generate_a2()
        hecke = h.HeckeAlgebra(group)

        left, right = hecke.generate_orders()

        # e, s, r, sr, rs, rsr
        expected_right = np.array([
            [1, 0, 0, 0, 0, 0],
            [1, 1, 0, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [1, 1, 0, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [1, 1, 1, 1, 1, 1]
        ], dtype=bool)
        expected_left = np.array([
            [1, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 1, 0],
            [1, 0, 1, 1, 0, 0],
            [1, 0, 1, 1, 0, 0],
            [1, 1, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1]
        ], dtype=bool)

        np.testing.assert_array_equal(right, expected_right)
        # np.testing.assert_array_equal(left, expected_left)
