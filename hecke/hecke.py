import hecke.laurent as l

class Hecke:
    def __init__(self, group, elements):
        """
        Initializes a hecke algebra.
        :param group: The Coxeter group.
        :param elements: A map from name to Laurent polynomial.
        """
        self.group = group
        self.elements = elements

    def __add__(self, other):
        ret = {}
        allkeys = set(self.elements.keys()).union(other.elements.keys())
        for key in allkeys:
            ret[key] = self.elements.get(key, l.zero) + other.elements.get(key, l.zero)

        return Hecke(self.group, ret)

    def __mul__(self, other):
        ret = Hecke(self.group, {})
        for thiselement, thiscoeff in self.elements.items():
            for otherelement, othercoeff in other.elements.items():
                ret += self.simple_mul(thiselement, thiscoeff,
                                       otherelement, othercoeff)
        return ret

    def simple_mul(self, thiselement, thiscoeff, otherelement, othercoeff):
        if otherelement == 'e' or otherelement == '':
            return Hecke(self.group, {thiselement: thiscoeff * othercoeff})
        else:
            s = otherelement[0]
            prod = self.simple_simple_mul(thiselement, thiscoeff, s, othercoeff)
            return prod * Hecke(self.group, {otherelement[1:]: l.one})

    def simple_simple_mul(self, thiselement, thiscoeff, s, othercoeff):
        w = self.group.get(thiselement)
        s = self.group.get(s)
        ws = w * s
        if ws.length() > w.length():
            return Hecke(self.group, {ws.name: thiscoeff * othercoeff})
        else:
            return Hecke(self.group, {
                ws.name: thiscoeff * othercoeff,
                w.name: thiscoeff * othercoeff * l.Laurent({-1: 1, 1: -1})
            })

    def __eq__(self, other):
        return self.group == other.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.group}, {self.elements})'
