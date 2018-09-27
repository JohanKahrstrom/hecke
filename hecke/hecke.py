import hecke.laurent as l
import copy


class HeckeAlgebra:
    def __init__(self, group):
        self.group = group
        self.zero = HeckeElement(self, {})
        self.one = HeckeElement(self, {'e': l.one})
        self._kl_basis = None
        self._dual_kl_basis = None

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
            for name in sorted(self.group.elements.keys(),
                               key=lambda str: (len(str), str)):
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

            self._kl_basis = ret
        return self._kl_basis

    def generate_dual_kl_basis(self):
        if self._dual_kl_basis is None:
            kl_basis = self.generate_kl_basis()
            ret = dict()

            # Add longest element
            ret[self.group.longest.name] = self.get_standard_basis_element(self.group.longest.name)

            # Add the rest
            for name in sorted(self.group.elements.keys(),
                               key=lambda str: (len(str), str),
                               reverse=True):
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

            self._dual_kl_basis = ret
        return self._dual_kl_basis

    def in_kl_basis(self, h):
        tmp = h.deepcopy()
        kl = self.generate_kl_basis()
        ret = dict()

        while tmp != self.zero:
            break

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

    def __getitem__(self, item):
        return self.elements.get(item, l.zero)

    def __eq__(self, other):
        return self.hecke.group == other.hecke.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.hecke.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.hecke.group}, {self.elements})'
