import sortedcontainers as sc
import collections


class Laurent:
    def __init__(self, coef):
        a = [(key, value) for key, value in coef.items() if value != 0]
        self.coef = sc.SortedDict(a)

    def __eq__(self, other):
        return self.coef == other.coef

    def __add__(self, other):
        c = collections.Counter()
        for deg1, value1 in self.coef.items():
            c[deg1] += value1
        for deg2, value2 in other.coef.items():
            c[deg2] += value2
        return Laurent(c)

    def __neg__(self):
        d = dict()
        for deg, value in self.coef.items():
            d[deg] = -value
        return Laurent(d)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, int):
            d = dict()
            for deg, value in self.coef.items():
                d[deg] = value * other
            return Laurent(d)
        else:
            c = collections.Counter()
            for deg1, value1 in self.coef.items():
                for deg2, value2 in other.coef.items():
                    c[deg1 + deg2] += value1 * value2
            return Laurent(c)

    def is_zero(self):
        return self == zero

    def involute(self):
        d = dict()
        for deg, value in self.coef.items():
            d[-deg] = value
        return Laurent(d)

    def __str__(self):
        return str(self.coef)

    def __repr__(self):
        return f'hecke.laurent.Laurent({self.coef})'


zero = Laurent({})
one = Laurent({0: 1})
