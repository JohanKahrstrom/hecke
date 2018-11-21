import sortedcontainers as sc
import collections


class Laurent:
    def __init__(self, coef):
        a = [(key, value) for key, value in coef.items() if value != 0]
        self.coef = sc.SortedDict(a)

    def shift(self, s):
        return Laurent({(key + s): value
                        for key, value in self.coef.items()})

    def copy(self):
        return Laurent(self.coef.copy())

    def __contains__(self, key):
        return key in self.coef

    def __getitem__(self, key):
        return self.coef.get(key, 0)

    def __call__(self, v):
        sum = 0
        for deg, coef in self.coef.items():
            sum += coef * (v ** deg)
        return sum

    def __iter__(self):
        for key, value in self.coef.items():
            yield (key, value)

    def items(self):
        return self.coef.items()

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

    def all_positive_degree(self):
        for deg, value in self.coef.items():
            if deg <= 0:
                return False
        return True

    def involute(self):
        d = dict()
        for deg, value in self.coef.items():
            d[-deg] = value
        return Laurent(d)

    def top(self):
        if len(self.coef) == 0:
            return None
        else:
            return self.coef.keys()[-1]

    def bottom(self):
        if len(self.coef) == 0:
            return None
        else:
            return self.coef.keys()[0]

    def shift(self, n):
        return self * Laurent({n: 1})

    def __str__(self):
        return str(self.coef)

    def __repr__(self):
        return f'hecke.laurent.Laurent({self.coef})'


zero = Laurent({})
one = Laurent({0: 1})
