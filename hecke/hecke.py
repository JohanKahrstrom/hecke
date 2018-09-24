import hecke.laurent as l


class HeckeAlgebra:
    def __init__(self, group):
        self.group = group

    def get_standard_basis_element(self, element):
        """Returns a standard basis element H_x"""
        if element not in self.group.elements:
            raise Exception(f"Can't create standard basis element for {element}.")
        return HeckeElement(self.group, {element: l.Laurent({0: 1})})


class HeckeElement:
    def __init__(self, hecke, elements):
        """
        Initializes a hecke algebra element.
        :param group: The Hecke algebra.
        :param elements: A map from name to Laurent polynomial.
        """
        self.hecke = hecke
        self.elements = elements

    def __add__(self, other):
        ret = {}
        allkeys = set(self.elements.keys()).union(other.elements.keys())
        for key in allkeys:
            ret[key] = self.elements.get(key, l.zero) + other.elements.get(key, l.zero)

        return HeckeElement(self.hecke, ret)

    def __mul__(self, other):
        ret = HeckeElement(self.hecke, {})
        for thiselement, thiscoeff in self.elements.items():
            for otherelement, othercoeff in other.elements.items():
                ret += self.simple_mul(thiselement, thiscoeff,
                                       otherelement, othercoeff)
        return ret

    def simple_mul(self, thiselement, thiscoeff, otherelement, othercoeff):
        if otherelement == 'e' or otherelement == '':
            return HeckeElement(self.hecke, {thiselement: thiscoeff * othercoeff})
        else:
            s = otherelement[0]
            prod = self.simple_simple_mul(thiselement, thiscoeff, s, othercoeff)
            return prod * HeckeElement(self.hecke, {otherelement[1:]: l.one})

    def simple_simple_mul(self, thiselement, thiscoeff, s, othercoeff):
        w = self.hecke.group.get(thiselement)
        s = self.hecke.group.get(s)
        ws = w * s
        if ws.length() > w.length():
            return HeckeElement(self.hecke, {ws.name: thiscoeff * othercoeff})
        else:
            return HeckeElement(self.hecke, {
                ws.name: thiscoeff * othercoeff,
                w.name: thiscoeff * othercoeff * l.Laurent({-1: 1, 1: -1})
            })

    def __eq__(self, other):
        return self.hecke.group == other.hecke.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.hecke.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.hecke.group}, {self.elements})'
