import hecke.laurent as l


class HeckeAlgebra:
    def __init__(self, group):
        self.group = group
        self.zero = HeckeElement(self, {})
        self.one = HeckeElement(self, {'e': l.one})

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
        w = self.group.get(thiselement)
        s = self.group.get(s)
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

    def __add__(self, other):
        ret = {}
        allkeys = set(self.elements.keys()).union(other.elements.keys())
        for key in allkeys:
            ret[key] = self.elements.get(key, l.zero) + other.elements.get(key, l.zero)

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

    def __eq__(self, other):
        return self.hecke.group == other.hecke.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.hecke.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.hecke.group}, {self.elements})'
