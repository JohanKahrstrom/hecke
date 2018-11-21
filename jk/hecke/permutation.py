import math

class Permutation:
    """
    values: A list encoding the permutation
    """
    def __init__(self, values):
        self.values = values

    def __eq__(self, other):
        return self.values == other.values

    def __mul__(self, other):
        def sgn(x):
            if x > 0:
                return 1
            elif x < 0:
                return -1
            else:
                return 0
        """
        'Composition on the right', i.e.
            p.Permutation([2, 1, 3]) *
            p.Permutation([1, 3, 2]) =
            p.Permutation([2, 3, 1])

        :param other:
        :return:
        """
        return Permutation([self.values[abs(v) - 1] * sgn(v) for v in other.values])

    def inverse(self):
        ret = self.values.copy()
        for index, element in enumerate(self.values):
            ret[element - 1] = index + 1
        return Permutation(ret)

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return f'hecke.permutation.Permutation({self.values})'

    def __str__(self):
        return str(self.values)

    def __hash__(self):
        hash_code = 1
        for e in self.values:
            hash_code = hash_code * 31 + hash(e)
        return hash_code
