"""

"""
import sys
sys.path.append('../')
import multiprocessing as ml
import time
import numpy as np
from linckage import msprime_simulate_variants
import collections as cl

# class Variant:
#     def __init__(self, genotype):
#         self.genotypes = genotype


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
    not set_genotype_1.intersection(set_genotype_2) \
    or set_genotype_1.union(set_genotype_2) == set_genotype_1 \
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
    list_incompatible = sorted(list_incompatible, key=lambda t: (t[0], t[1])) #
    index = len(list_incompatible) - 1
    while index > 0:
        min1, max1 = list_incompatible[index - 1]
        min2, max2 = list_incompatible[index]
        if min1 <= min2 and max1 >= max2: #if the second tuple is emcompassed within the first
            list_incompatible[index - 1] = list_incompatible[index]
            del list_incompatible[index]
        elif min1 <= min2 and max1 <= max2 and max1 > min2: # if the two tuples are overlapping
            del list_incompatible[index]
            list_incompatible[index - 1] = (min2, max1)
        index -= 1
    return list_incompatible


def detect_internal_incompatibilities(variants, oriented=True, thresold=10):
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


def accuracy(ts, th, unique=False, oriented=True):
    variants = [variant for variant in ts.variants() if sum(variant.genotypes) > 1]
    breakpoints = list(ts.breakpoints())
    mld = detect_internal_incompatibilities(variants, oriented, th)
    nb_recombinaison = len(breakpoints)
    nb_observed = len(mld)
    count = 0
    for breakpoint in breakpoints:
        for min, max in mld:
            if breakpoint > variants[min].site.position and \
            breakpoint < variants[max].site.position:
                count += 1
                if unique:
                    mld.remove((min, max))
                break
    return count / nb_recombinaison if not unique \
    else (count / nb_recombinaison, len(mld) / nb_observed)

def breakpoint_by_incompatility(ts, th, oriented=True):
    variants = [variant for variant in ts.variants() if sum(variant.genotypes) > 1]
    breakpoints = list(ts.breakpoints())
    mld = detect_internal_incompatibilities(variants, oriented, th)
    dict_break_mld = {key:0 for key in mld}
    for breakpoint in breakpoints:
        for min, max in mld:
            if breakpoint > variants[min].site.position and \
            breakpoint < variants[max].site.position:
                dict_break_mld[(min, max)] += 1
               # mld.remove((min, max))
                break
    return cl.Counter(dict_break_mld.values())


def main():
    start = time.time()
    params = {"sample_size": 20, "Ne": 1, "ro": 8e-4, "mu": 8e-3,  "Tau": 1.0,
    "Kappa": 1.0 , "length": int(1e5)}
    ts = msprime_simulate_variants(params)
    accuracy(ts, True)
    accuracy(ts, oriented=False)


if __name__ == "__main__":
    main()
