import time
import numpy as np
import itertools as it
from msprime_simulation import msprime_simulate_variants

def internal_incompatibility(genotype1, genotype2, oriented=True):
    """
    four gamete test enlarged to discrete recombnation events since genotypes
    are oriented ( 0 = ancestral state, 1 = derivative state)
    two genotypes are compatible, there is no breakpoint between them,
    when they exlude each others or one of them include the other
    parameters:
        genotype1, ganotype2, array of int values taken only 0 or 1 as value
    return:
        boolean, True if the two sites if there is an incompatible or discrete
        recombination event between them, then False
    """
    if not oriented:
        return len(set([*zip(genotype1, genotype2)])) == 4
    set_genotype_1 = {i for i in range(len(genotype1)) if genotype1[i]}
    set_genotype_2 = {i for i in range(len(genotype2)) if genotype2[i]}
    return not (
        not set_genotype_1.intersection(set_genotype_2)
        or set_genotype_1.union(set_genotype_2) == set_genotype_1
        or set_genotype_1.union(set_genotype_2) == set_genotype_2
    )


def shortening(list_incompatible):
    """
    parameters:
    list_incompatible, list of int's tuples, each tuple define two incompatible variant
    sites defining a breakpoint
    return:
    list_incompatible, list of non overlapping int's tuples
    """
    list_incompatible = sorted(list_incompatible, key=lambda t: (t[0], t[1]))
    index = len(list_incompatible) - 1
    while index > 0:
        min1, max1 = list_incompatible[index - 1]
        min2, max2 = list_incompatible[index]
        if min1 <= min2 and max1 >= max2:  # if the second tuple is emcompassed within the first
            list_incompatible[index - 1] = list_incompatible[index]
            del list_incompatible[index]
        elif (min1 <= min2 < max1) and max1 < max2:  # if the two tuples are overlapping
            del list_incompatible[index]
            list_incompatible[index - 1] = (min2, max1)
        index -= 1
    return list_incompatible


def choping(mld, variants, length):
    """
    entrée:
        mld: liste de tuples définissant un point de recombinaisons
        variants: liste d'objet Variants du module msprime
        length: longueur de l'alignement
    sortie:
        liste de float: position des points de recombinaisons sur l'alignement
    """
    return [0.0] + [(variants[end].site.position + variants[start].site.position) / 2
                    for start, end in mld] + [length]


def detect_internal_incompatibilities(variants, oriented=True, thresold=20):
    """
    detection of the incomaptible paires of variant sites between two born considering
    only sequences displaying derivative state on the segregative site.
    paramters:
        variants, list of msprime Variant class's instances,
        segregated, boolean list, true if, for the segragarive site, the corresponding
        sequence display the derivative state.
        born_inf, int: the closest site incompatible site with the segregative_site
        before him
        born_sup, int: the closest site incompatible site with the segregative_site
        after
        thresold, int:
    """
    list_incompatible_sites = []
    for i in range(len(variants)):
        j = i + 1
        while j < i + thresold and j < len(variants):
            if internal_incompatibility(variants[i].genotypes,
               variants[j].genotypes, oriented):
                list_incompatible_sites.append((i, j))
            j += 1
    return shortening(list_incompatible_sites)


def internal_incompatibility_2(set_genotype_1, set_genotype_2):
    """
    four gamete test enlarged to discrete recombnation events since genotypes
    are oriented ( 0 = ancestral state, 1 = derivative state)
    two genotypes are compatible, there is no breakpoint between them,
    when they exlude each others or one of them include the other
    parameters:
        genotype1, ganotype2, array of int values taken only 0 or 1 as value
    return:
        boolean, True if the two sites if there is an incompatible or discrete
        recombination event between them, then False
    """
    return not (
        not set_genotype_1.intersection(set_genotype_2)
        or set_genotype_1.union(set_genotype_2) == set_genotype_1
        or set_genotype_1.union(set_genotype_2) == set_genotype_2
    )


def closest_incompatibility(index: object, sub_set: object) -> object:
    position = [-1]
    for sub_set_tmp in index:
        if internal_incompatibility_2(set(sub_set), set(sub_set_tmp)):
            position.append(index[sub_set_tmp])
    return max(position)


def built_index(start, max_start, variants, individuals):
    block_start, position = start, start
    index = {tuple(individuals[variants[start].genotypes == 1]): start}
    while position < len(variants):
        variant = variants[position]
        sub_set = tuple(individuals[variant.genotypes == 1])
        if len(sub_set) != 1:
            if sub_set not in index:
                position_2 = closest_incompatibility(index, sub_set)
                if position_2 != -1:
                    if position_2 > max_start:
                        return block_start, position_2, position
                    block_start, position = position_2 + 1, position_2 + 1
                    sub_set = tuple(individuals[variants[position].genotypes == 1])
                    index = {sub_set: position}
                    position += 1
                    continue
            index[sub_set] = position
        position += 1
    return block_start, - 1, len(variants) - 1


def detect_events(variants, nb):
    individuals = np.array((range(nb)))
    block_start, start, block_end = built_index(0, 0, variants, individuals)
    list_block = [(block_start, block_end)]
    while block_end < len(variants) - 1:
        block_start, start, block_end = built_index(start + 1, block_end, variants, individuals)
        list_block.append((block_start - 1, block_end))
    return [(list_block[i + 1][0], list_block[i][1]) for i in range(len(list_block) - 1)] #list_block


# def detect_events(variants, nb):
#     index = {}
#     list_mld = []
#     individuals = np.array((range(nb)))
#     for position, variant in enumerate(variants):
#         sub_set = tuple(individuals[variant.genotypes == 1])
#         if len(sub_set) != 1:
#             if sub_set not in index:
#                 position_2 = closest_incompatibility(index, sub_set)
#                 if position_2 != -1:
#                     print(index, len(index))
#                     index = {sub_set: position}
#                     list_mld.append((position_2, position))
#                     continue
#             index[sub_set] = position
#     return list_mld


def main():
    params = {"sample_size": 10, "Ne": 1, "ro": 8e-6, "mu": 8e-4, "Tau": 1, "Kappa": 1, "length": int(1e7)}
    _, _, _, variants = msprime_simulate_variants(params)
    #print(detect_events(variants, 10))
    # print(len(detect_events(variants, 10)))
    # print(len(detect_internal_incompatibilities(variants, oriented=True, thresold=20)))
    # print(time.time() - start)

main()