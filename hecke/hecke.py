import hecke.laurent as l
import copy
import numpy as np
import time


class HeckeAlgebra:
    def __init__(self, group):
        self.group = group
        self.zero = HeckeElement(self, {})
        self.one = HeckeElement(self, {'e': l.one})
        self._kl_basis = None
        self._dual_kl_basis = None
        self._left_order = None
        self._right_order = None

    def element(self, d):
        """
        Creates a hecke element given by the map d.
        :param d:
        :return:
        """
        return HeckeElement(self, d)

    def simple_mul(self, thiselement, otherelement):
        """
        Calculates the product H_x*H_y.
        """
        if otherelement == 'e' or otherelement == '':
            return self.get_standard_basis_element(thiselement)
        else:
            s = otherelement[0]
            prod = self.simple_simple_mul(thiselement, s)
            return prod * HeckeElement(self, {otherelement[1:]: l.one})

    def simple_simple_mul(self, thiselement, s):
        w = self.group[thiselement]
        s = self.group[s]
        ws = w * s
        if ws.length() > w.length():
            return HeckeElement(self, {ws.name: l.one})
        else:
            return HeckeElement(self, {
                ws.name: l.one,
                w.name: l.Laurent({-1: 1, 1: -1})
            })

    def get_standard_basis_element(self, element):
        """Returns a standard basis element H_x"""
        if element not in self.group.elements:
            raise Exception(f"Can't create standard basis element for {element}.")
        return HeckeElement(self, {element: l.one})

    def get_standard_inverse_element(self, element):
        if element not in self.group.elements:
            raise Exception(f"Can't create standard inverse element for {element}.")
        ret = self.one
        for generator in element[::-1]:
            ret *= self.get_generator_inverse_element(generator)
        return ret

    def get_generator_inverse_element(self, generator):
        return self.element({
            'e': l.Laurent({-1: -1, 1: 1}),
            generator: l.one
        })

    def get_standard_dual(self, element):
        if element not in self.group.elements:
            raise Exception(f"Can't create standard inverse element for {element}.")
        ret = self.one
        for generator in element:
            ret *= self.get_generator_inverse_element(generator)
        return ret

    def generate_kl_basis(self):
        if self._kl_basis is None:
            n = len(self.group.all_elements())
            ret = dict()

            # Add identity
            ret['e'] = self.get_standard_basis_element('e')

            # Add generators
            for generator in self.group.generators:
                ret[generator] = (
                        self.get_standard_basis_element(generator) +
                        self.get_standard_basis_element('e') * l.Laurent({1: 1})
                )

            # Add the rest
            start = time.time()
            for i, name in enumerate(sorted(self.group.elements.keys(),
                                            key=lambda str: (len(str), str))):
                if name not in ret:
                    x = name[:-1]
                    kl_x = ret[x]
                    s = name[-1]
                    kl_s = ret[s]
                    # Now x < xs
                    xs = name
                    kl_xs = kl_x * kl_s
                    sub = self.zero
                    if xs not in ret:
                        # Subtract all y[1] with ys < y:
                        for y, coef in kl_x.elements.items():
                            if 1 in coef and (self.group[y] * self.group[s]).length() < self.group[y].length():
                                sub -= ret[y] * coef[1]
                        ret[xs] = kl_xs + sub
                    difftime = time.time() - start
                    if difftime > 5.0:
                        print('KL-basis generated ' + '{:2.2f}%'.format(
                            100 * (i + 1) / n))
                        start = time.time()

            self._kl_basis = ret
        return self._kl_basis

    def generate_dual_kl_basis(self):
        if self._dual_kl_basis is None:
            kl_basis = self.generate_kl_basis()
            n = len(self.group.all_elements())
            ret = dict()

            # Add longest element
            ret[self.group.longest.name] = self.get_standard_basis_element(self.group.longest.name)

            start = time.time()
            # Add the rest
            for i, name in enumerate(sorted(self.group.elements.keys(),
                                            key=lambda str: (len(str), str),
                                            reverse=True)):
                x = self.group[name]
                for generator in self.group.generators.keys():
                    s = self.group[generator]
                    xs = x * s
                    if xs.length() < x.length() and xs.name not in ret:
                        dkl_x = ret[x.name]
                        kl_s = kl_basis[generator]
                        dkl_xs = dkl_x * kl_s
                        sub = dkl_x * l.Laurent({-1: 1, 1: 1})
                        # Subtract all y[1] with ys > y:
                        for y, coef in dkl_x.elements.items():
                            if 1 in coef and (self.group[y] * self.group[s.name]).length() > self.group[y].length():
                                sub -= ret[y] * coef[1]
                        ret[xs.name] = dkl_xs - sub
                difftime = time.time() - start
                if difftime > 5.0:
                    print('Dual KL-basis generated ' + '{:2.2f}%'.format(100 * (i + 1) / n))
                    start = time.time()

            self._dual_kl_basis = ret
        return self._dual_kl_basis

    def generate_orders(self):
        if self._left_order is None or self._right_order is None:
            digraph = self._generate_digraph()
            self._right_order = self._order_from_digraph(digraph)

            # There's probably a better way of doing this...
            self._left_order = np.zeros(self._right_order.shape)
            for xi, x in enumerate(self.group.all_elements()):
                xinversei = x.inverse().index
                for yi, y in enumerate(self.group.all_elements()):
                    yinversei = y.inverse().index
                    self._left_order[xi, yi] = self._right_order[xinversei, yinversei]

        return self._left_order, self._right_order

    def _generate_digraph(self):
        """
        Returns a map element -> set of elements directly larger.
        :return:
        """
        kl_basis = self.generate_kl_basis()
        print('Generating digraph')
        n = len(self.group.all_elements())
        d = dict()
        start = time.time()
        for i, x in enumerate(self.group.all_elements()):
            s = set()
            for g in self.group.all_generators():
                s.update((kl_basis[x.name] * kl_basis[g.name]).in_kl_basis().keys())
            d[x.name] = s
            difftime = time.time() - start
            if difftime > 5.0:
                print('Digram generated ' + '{:2.2f}%'.format(
                    100 * (i + 1) / n))
                start = time.time()
        print('Finished generating digraph')
        return d

    def _order_from_digraph(self, digraph):
        def rec(base, x, order):
            for yname in digraph[x.name]:
                y = self.group[yname]
                if not order[y.index, base.index]:
                    order[y.index, base.index] = True
                    rec(base, y, order)

        print('Generating order from digraph')
        n = len(self.group.elements)
        order = np.diag(np.ones(n, dtype=bool))
        order[self.group.identity.index, self.group.identity.index] = True
        start = time.time()
        for i, x in enumerate(self.group.all_elements()[::-1]):
            rec(x, x, order)
            difftime = time.time() - start
            if difftime > 5.0:
                print('Order from digraph generated ' + '{:2.2f}%'.format(
                    100 * (i + 1) / n))
                start = time.time()
        print('Finished generating order from digraph')
        return order

    def filtration_str(self, d):
        """Returns a string repreentation of the filtration 'd'. 'd' is a
        dictionary from element name to Laurent polynomial"""
        def get_coefficient_string(value):
            if value == 1:
                return ""
            else:
                return str(value)

        if not d:
            return "\n"

        by_degree = dict()
        for name in sorted(d.keys(), key=lambda x: (len(x), x)):
            coef = d[name]
            for degree, value in coef.items():
                if degree in by_degree:
                    by_degree[degree] += f" {get_coefficient_string(value)}{name}"
                else:
                    by_degree[degree] = f"{get_coefficient_string(value)}{name}"

        if not by_degree:
            return "\n"

        width = max([len(s) for s in by_degree.values()])
        bottom = min(by_degree.keys())
        top = max(by_degree.keys())
        return "\n".join([
            by_degree.get(deg, "").center(width, " ")
            for deg in range(bottom, top + 1)
        ]) + "\n"

    def basis_matrix(self, elements, basis):
        """
        Returns a matrix representing the change of basis between basis1 and
        basis2, when specializing v = 1.
        :param basis1:
        :param basis2:
        :return:
        """
        d = dict()
        for name, x in elements:
            row = dict()
            xb = x.in_basis(basis)
            for y in sorted(basis.keys(), key=lambda x: (len(x), x)):
                p = xb.get(y, l.zero)
                row[y] = p
            d[name] = row
        return self.matrix(d)

    def vector(self, d):
        """d is a dictionary element -> Laurent"""
        ret = []
        for x in self.group.all_elements():
            p = d.get(x.name, l.zero)
            ret.append(p(1))
        return ret

    def matrix(self, d):
        """d is a dictionary Dict[element, Dict[element, value]]"""
        ret = []
        for x in self.group.all_elements():
            row = []
            xr = d.get(x.name, dict())
            ret.append(self.vector(xr))
        return ret


class HeckeElement:
    def __init__(self, hecke, elements):
        """
        Initializes a hecke algebra element.
        :param group: The Hecke algebra.
        :param elements: A map from name to Laurent polynomial.
        """
        self.hecke = hecke
        self.elements = {element: coeff for element, coeff in elements.items()
                         if coeff != l.zero}

    def deepcopy(self):
        return HeckeElement(self.hecke, copy.deepcopy(self.elements))

    def __add__(self, other):
        ret = {}
        allkeys = set(self.elements.keys()).union(other.elements.keys())
        for key in allkeys:
            ret[key] = self.elements.get(key, l.zero) + other.elements.get(key, l.zero)

        return HeckeElement(self.hecke, ret)

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        ret = {}
        for element, coeff in self.elements.items():
            ret[element] = -coeff
        return HeckeElement(self.hecke, ret)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, l.Laurent):
            return self.multiply_int(other)
        else:
            return self.multiply_hecke(other)

    def multiply_int(self, other):
        """Also works if other is l.Laurent"""
        ret = self.hecke.zero
        for thiselement, thiscoeff in self.elements.items():
            ret += HeckeElement(self, {thiselement: thiscoeff * other})
        return ret

    def multiply_hecke(self, other):
        ret = self.hecke.zero
        for thiselement, thiscoeff in self.elements.items():
            for otherelement, othercoeff in other.elements.items():
                coeff = thiscoeff * othercoeff
                ret += self.hecke.simple_mul(thiselement, otherelement) * coeff
        return ret

    def shift(self, n):
        return self * l.Laurent({n: 1})

    def dual(self):
        ret = self.hecke.zero
        for element, coef in self.elements.items():
            ret += self.hecke.get_standard_dual(element) * coef.involute()
        return ret

    def i(self):
        """H_x.i() = H_x^-1"""
        ret = self.hecke.zero
        for element, coef in self.elements.items():
            ret += (
                self.hecke.get_standard_basis_element(
                    self.hecke.group[element].inverse().name
                ) * coef
            )
        return ret

    def tau(self):
        return self['e']

    def bottom(self):
        """Returns the set of keys that have lowest degree coefficient"""
        ret = []
        mindegree = None
        for element, coef in self.elements.items():
            if mindegree == None or mindegree > coef.bottom():
                mindegree = coef.bottom()
                ret = [(element, mindegree)]
            elif mindegree == coef.bottom():
                ret.append((element, mindegree))

        return ret

    def __getitem__(self, item):
        return self.elements.get(item, l.zero)

    def __eq__(self, other):
        return self.hecke.group == other.hecke.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.hecke.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.hecke.group}, {self.elements})'

    def in_basis(self, basis):
        """
        Decomposes 'h' into the basis 'basis'. 'basis' is a dictionary from
        element name to Laurent polynomial, which all have to be
        :param basis:
        :param h:
        :return:
        """
        tmp = self.deepcopy()
        ret = dict()

        while tmp != self.hecke.zero:
            for bottom_element, degree in tmp.bottom():
                coeff = tmp[bottom_element][degree]
                sub = basis[bottom_element].shift(degree) * coeff
                ret[bottom_element] = (
                        l.Laurent({degree: coeff}) +
                        ret.get(bottom_element, l.zero)
                )
                tmp -= sub

        return ret

    def in_kl_basis(self):
        basis = self.hecke.generate_kl_basis()
        return self.in_basis(basis)

    def in_dual_kl_basis(self):
        basis = self.hecke.generate_dual_kl_basis()
        return self.in_basis(basis)

    def dual_kl_filtration(self):
        return self.hecke.filtration_str(self.in_dual_kl_basis())
