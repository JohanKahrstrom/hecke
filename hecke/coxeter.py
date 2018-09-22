import hecke.permutation as p
import sortedcontainers as sc

"""
Perhaps not the ideal set-up. We have two linked classes, where the element
only knows its name and which group it belongs to. The group knows which
permutations each name corresponds to, as well as a map from permutation
to name.
"""
class CoxeterElement:
    def __init__(self, group, name):
        self.group = group
        self.name = name

    def __mul__(self, other):
        if self.group != other.group:
            raise Exception("Can't multiply elements from different groups")
        # Calculate the new permutation
        new_permutation = self.group.permutations[self.name] * other.group.permutations[other.name]
        # Find the name of this element
        new_name = self.group.name_lookup[new_permutation]
        return self.group.get(new_name)

    def length(self):
        if self.name == 'e':
            return 0
        else:
            return len(self.name)

    def __eq__(self, other):
        return self.group == other.group and self.name == other.name


class CoxeterGroup:
    def __init__(self, elements, permutations, name_lookup):
        self.elements = elements
        self.permutations = permutations
        self.name_lookup = name_lookup

    def get(self, name):
        return self.elements[name]

    def get_x(self, name):
        element = self.get('e')
        for generator_name in list(name):
            element *= self.get(generator_name)
        return element

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
        group = CoxeterGroup(elements,
                             permutations,
                             name_lookup)

        lengths = set([len(v) for v in generators.values()])
        if len(lengths) != 1:
            raise Exception(f'Generators should all have the same length (lengths: {lengths})')

        seen_permutations = set()
        # Add identity
        identity_permutation = p.Permutation(list(range(1, list(lengths)[0] + 1)))
        elements['e'] = CoxeterElement(group, 'e')
        group.identity = elements['e']
        permutations['e'] = identity_permutation
        seen_permutations.add(identity_permutation)
        # Add generators
        for name, permutation in sorted_generators.items():
            elements[name] = CoxeterElement(group, name)
            seen_permutations.add(permutation)
            permutations[name] = permutation
        # Generate the rest of the elements
        last_layer = generators
        while len(last_layer) > 0 and len(seen_permutations) < 10:
            new_permutations = dict()
            for name, permutation in last_layer.items():
                for gname, generator in sorted_generators.items():
                    new_name = name + gname
                    new_permutation = permutation * generator
                    if not new_permutation in seen_permutations:
                        permutations[new_name] = new_permutation
                        seen_permutations.add(new_permutation)
                        elements[new_name] = CoxeterElement(group, new_name)
                        new_permutations[new_name] = new_permutation
            last_layer = new_permutations

        for name, permutation in permutations.items():
            name_lookup[permutation] = name

        return group


def generate_a1():
    return CoxeterGroup.generate({'s': p.Permutation([2, 1])})


def generate_a2():
    return CoxeterGroup.generate({
            'r': p.Permutation([2, 1, 3]),
            's': p.Permutation([1, 3, 2])
        })