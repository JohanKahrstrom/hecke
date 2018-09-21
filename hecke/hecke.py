import collections
import hecke.laurent as l

class Hecke:
    def __init__(self, group, elements):
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
                if otherelement == 'e':
                    ret += Hecke(self.group, {thiselement: thiscoeff * othercoeff})
                else
                    pass

        return Hecke(self.group, ret)

    def __eq__(self, other):
        return self.group == other.group and self.elements == other.elements

    def __repr__(self):
        return f'hecke.hecke.Hecke({self.group}, {self.elements})'

    def __str__(self):
        return f'Hecke({self.group}, {self.elements})'
