import hecke.permutation as p
import sortedcontainers as sc

"""
Perhaps not the ideal set-up. We have two linked classes, where the element
only knows its name and which group it belongs to. The group knows which
permutations each name corresponds to, as well as a map from permutation
to name.
"""
class CoxeterElement:
    def __init__(self, group, name, index):
        self.group = group
        self.name = name
        self.index = index

    def __mul__(self, other):
        if self.group != other.group:
            raise Exception("Can't multiply elements from different groups")
        # Calculate the new permutation
        new_permutation = self.group.permutations[self.name] * other.group.permutations[other.name]
        # Find the name of this element
        new_name = self.group.name_lookup[new_permutation]
        return self.group[new_name]

    def inverse(self):
        return self.group.get_inverse(self.name)

    def permutation(self):
        return self.group.permutations[self.name]

    def length(self):
        if self.name == 'e':
            return 0
        else:
            return len(self.name)

    def __eq__(self, other):
        return self.group == other.group and self.name == other.name


class CoxeterGroup:
    def __init__(self, generators, elements, permutations, name_lookup, inverse):
        self.generators = generators
        self._all_generators = None
        self.elements = elements
        self._all_elements = None
        self.permutations = permutations
        self.name_lookup = name_lookup
        self.inverse = inverse

    def __getitem__(self, name):
        return self.elements[name]

    def get_x(self, name):
        element = self['e']
        for generator_name in list(name):
            element *= self.get(generator_name)
        return element

    def get_inverse(self, name):
        return self.inverse[name]

    @staticmethod
    def generate(generators):
        """
        Generates a group from a dictionary of generators
        { name: permutation }.
        :param generators:
        :return:
        """
        sorted_generators = sc.SortedDict(generators)
        elements = dict()
        permutations = dict()
        name_lookup = dict()
        inverses = dict()
        group = CoxeterGroup(generators,
                             elements,
                             permutations,
                             name_lookup,
                             inverses)

        lengths = set([len(v) for v in generators.values()])
        if len(lengths) != 1:
            raise Exception(f'Generators should all have the same length (lengths: {lengths})')

        seen_permutations = set()
        # Add identity
        identity_permutation = p.Permutation(list(range(1, list(lengths)[0] + 1)))
        elements['e'] = CoxeterElement(group, 'e', 0)
        group.identity = elements['e']
        permutations['e'] = identity_permutation
        seen_permutations.add(identity_permutation)
        # Add generators
        for name, permutation in sorted_generators.items():
            elements[name] = CoxeterElement(group, name, len(elements))
            seen_permutations.add(permutation)
            permutations[name] = permutation
        # Generate the rest of the elements
        last_layer = generators
        nrloops = 0
        while len(last_layer) > 0:
            new_permutations = dict()
            for name, permutation in last_layer.items():
                for gname, generator in sorted_generators.items():
                    new_name = name + gname
                    new_permutation = permutation * generator
                    if not new_permutation in seen_permutations:
                        permutations[new_name] = new_permutation
                        seen_permutations.add(new_permutation)
                        elements[new_name] = CoxeterElement(group, new_name, len(elements))
                        new_permutations[new_name] = new_permutation
            last_layer = new_permutations
            nrloops += 1

        for name, permutation in permutations.items():
            name_lookup[permutation] = name

        # Calculate inverses
        for name, permutation in permutations.items():
            inverse = identity_permutation
            inverse_name = name[::-1] # Reversed string
            for generator in inverse_name:
                inverse = inverse * permutations[generator]
            inverses[name] = elements[name_lookup[inverse]]

        # Add identity and longest element
        group.identity = elements['e']
        longest_name = sorted(elements.keys(), key=lambda x: (len(x), x))[-1]
        group.longest = elements[longest_name]

        return group

    def all_generators(self):
        """Returns the generators, ordered by name"""
        if self._all_generators is None:
            self._all_generators = [self.elements[key]
                                    for key
                                    in sorted(self.generators.keys())]
        return self._all_generators

    def all_elements(self):
        """Returns the elements, ordered by (len(nate), name)"""
        if self._all_elements is None:
            self._all_elements = [self.elements[key]
                                  for key
                                  in sorted(self.elements.keys(), key=lambda x: (len(x), x))]
        return self._all_elements


def generate_a1():
    return CoxeterGroup.generate({'s': p.Permutation([2, 1])})


def generate_a2():
    return CoxeterGroup.generate({
            'r': p.Permutation([2, 1, 3]),
            's': p.Permutation([1, 3, 2])
        })


def generate_a3():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4]),
        's': p.Permutation([1, 3, 2, 4]),
        't': p.Permutation([1, 2, 4, 3])
    })


def generate_a4():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4, 5]),
        's': p.Permutation([1, 3, 2, 4, 5]),
        't': p.Permutation([1, 2, 4, 3, 5]),
        'u': p.Permutation([1, 2, 3, 5, 4])
    })


def generate_a5():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4, 5, 6]),
        's': p.Permutation([1, 3, 2, 4, 5, 6]),
        't': p.Permutation([1, 2, 4, 3, 5, 6]),
        'u': p.Permutation([1, 2, 3, 5, 4, 6]),
        'v': p.Permutation([1, 2, 3, 4, 6, 5])
    })


def generate_b2():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1]),
        's': p.Permutation([1, -2])
    })


def generate_b3():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3]),
        's': p.Permutation([1, 3, 2]),
        't': p.Permutation([1, 2, -3])
    })


def generate_b4():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4]),
        's': p.Permutation([1, 3, 2, 4]),
        't': p.Permutation([1, 2, 4, 3]),
        'u': p.Permutation([1, 2, 3, -4])
    })


def generate_b5():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4, 5]),
        's': p.Permutation([1, 3, 2, 4, 5]),
        't': p.Permutation([1, 2, 4, 3, 5]),
        'u': p.Permutation([1, 2, 3, 5, 4]),
        'v': p.Permutation([1, 2, 3, 4, -5])
    })


def generate_d4():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4]),
        's': p.Permutation([1, 3, 2, 4]),
        't': p.Permutation([1, 2, 4, 3]),
        'u': p.Permutation([1, 2, -4, -3])
    })


def generate_d5():
    return CoxeterGroup.generate({
        'r': p.Permutation([2, 1, 3, 4, 5]),
        's': p.Permutation([1, 3, 2, 4, 5]),
        't': p.Permutation([1, 2, 4, 3, 5]),
        'u': p.Permutation([1, 2, 3, 5, 4]),
        'v': p.Permutation([1, 2, 3, -5, -4])
    })


def generate_g2():
    return CoxeterGroup.generate({
        'r': p.Permutation([1, 1]),
        's': p.Permutation([1, -2])
    })

# G_2
#      Coxeter.set_identity(CoxeterTypeGRepresentation.new([0,  1]), "e")
# r = Coxeter.add_generator(CoxeterTypeGRepresentation.new([0, -1]), "r")
# s = Coxeter.add_generator(CoxeterTypeGRepresentation.new([1, -1]), "s")