"""

"""
import multiprocessing as ml
import time
import numpy as np
from linckage import msprime_simulate_variants


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
    set_genotype_1 = {i for i in range(len(genotype1)) if genotype1[i]}
    set_genotype_2 = {i for i in range(len(genotype2)) if genotype2[i]}
    return not (
    not set_genotype_1.intersection(set_genotype_2) \
    or set_genotype_1.union(set_genotype_2) == set_genotype_1 \
    or set_genotype_1.union(set_genotype_2) == set_genotype_2
    )


def external_incompatibility(segregative_site, variants, revert=False):
    """
    test the presence of discrete or incompatible breakpoint between a segregative
    site and his neighbors
    parameters:
        segregative_site, tuple including a list of boolean,
            true if the coresponding sequence expose the derivative state,
            and an int, the positoin of the coresponding variant site in the list of
            variant.
        variants, list of msprime Variant class's instances
        revert, boolean: true if we look for an incompatible site toward the left in the
            variants list, false for the rigth
    return, int: the position of the closest incompatible site of the segregative site,
    if there is no incompatible site the return value is -1
    """
    segregated, position = segregative_site
    if revert:
        position = len(variants) - position
    for index, variant in enumerate(variants[position:]):
        # the current site are incompatible with the segregative site
        if len(set(variant.genotypes[segregated])) == 2 \
        and 1 in variant.genotypes[np.invert(segregated)]:
            return  len(variants[position:]) - index - 1 if revert else index + position
    # return  0 if revert else len(variants) - 1
    return -1


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
        # print(list_incompatible)
        min1, max1 = list_incompatible[index - 1]
        min2, max2 = list_incompatible[index]
        if min1 <= min2 and max1 >= max2: #if the second tuple is emcompassed within the first
            list_incompatible[index - 1] = list_incompatible[index]
            del list_incompatible[index]
        # elif min1 >= min2 and max1 <= max2:# #if the first tuple is emcompassed within the second
        #     del list_incompatible[index]
        elif min1 <= min2 and max1 <= max2 and max1 > min2: # if the two tuples are overlapping
            del list_incompatible[index]
            list_incompatible[index - 1] = (min2, max1)
        # elif min2 <= min1 and min2 >= max2 and max2 <= max1:
        #     del list_incompatible[index]
        #     list_incompatible[index - 1] = (min1, max2)
        index -= 1
    return list_incompatible


def segregative_sites(variants):
    """
    parameters:
    variants, list of msprime Variant class's instances
    return:
    list of tuple encommassing a tuple of boolean, true if the corresponding sequence deisplay
     with between three and sample size - 1 sequences displaying the derivative site
    """
    list_segregative_sites = []
    for index, variant in enumerate(variants):
        tmp = sum(variant.genotypes)
        if tmp > 2 and tmp < len(variant.genotypes) :
            list_segregative_sites.append((list(variant.genotypes == 1), index))
    return list_segregative_sites


def detect_internal_incompatibilities(variants, segregated, borne_inf, borne_sup, thresold=150):
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
    for i in range(borne_inf, borne_sup + 1):
        j = i + 1
        while j < i + thresold and j < borne_sup + 1:
            if internal_incompatibility(variants[i].genotypes[segregated],
            variants[j].genotypes[segregated]):
                list_incompatible_sites.append((i, j))
            j += 1
    return list_incompatible_sites


def non_overlaping_incompatibilities(variants, segregative_site):
    """
    paramters:
        variants, list of msprime Variant class's instances
        segregative_site, tuple including a list of boolean,
            true if the coresponding sequence expose the derivative state,
            and an int, the positoin of the coresponding variant site in the list of
            variant.
    return:
    list of non overlapping int's tuples
    """
    segregated, position_segragative_site = segregative_site
    #closest incomaptible site after, for borne_sup, and before, for borne_sup,
    #the segregative site
    borne_sup = external_incompatibility(segregative_site, variants, revert=False)
    borne_inf = external_incompatibility(segregative_site, variants[::-1], revert=True)
    #if there is no incompatible site before the segragative site
    borne_inf_tmp = 0 if borne_inf == -1 else borne_inf
    #if there is no incompatible site after the segragative site
    borne_sup_tmp = len(variants) - 1 if borne_sup == -1 else borne_sup
    #detection af all paires of incompatible sites between borne_sup and borne_inf
    list_incompatible_sites = detect_internal_incompatibilities(variants, segregated,
    borne_inf_tmp, borne_sup_tmp)
    if born_inf > 0:
        list_incompatible_sites.append((borne_inf, position_segragative_site))
    if born_sup > 0:
        list_incompatible_sites.append((position_segragative_site, borne_sup))
    return shortening(list_incompatible_sites)


def detect_mld_block(variants, segregated_site):
    """
    detect the mld block encompassing a given segregative site considering recombination
    events which implying at least one sequence displaying the derivative state at the considered
    site.
    paramters:
        variants, list of msprime Variant class's instances
        segregative_site, tuple including a list of boolean,
            true if the coresponding sequence expose the derivative state,
            and an int, the positoin of the coresponding variant site in the list of
            variant.
    """
    segregated, position_segragative_site = segregated_site
    list_incompatible_sites = non_overlaping_incompatibilities(variants, segregated_site)
    for index in range(len(list_incompatible_sites) - 1):
        borne_sup_inf = list_incompatible_sites[index][1]
        borne_inf_sup = list_incompatible_sites[index + 1][0]
        if borne_sup_inf <= position_segragative_site \
        and borne_inf_sup >= position_segragative_site:
            return (list_incompatible_sites[index], \
            list_incompatible_sites[index + 1], \
            segregated)
    return((0, 0), (len(variants) - 1, len(variants) -1), segregated)


def mld_blocks_for_all_segregative_sites(variants):
    list_segregative_sites = segregative_sites(variants)
    with ml.Pool(4) as pool:
        data = pool.starmap(detect_mld_block, [(variants, segregated) \
        for segregated in list_segregative_sites])
    return data


def main():
    start = time.time()
    params = {"sample_size": 8, "Ne": 1, "ro": 8e-4, "mu": 8e-2,  "Tau": 1.0,
    "Kappa": 1.0 , "length": int(1e5)}
    variants = list(msprime_simulate_variants(params).variants())
    print(len(variants))
    print(mld_blocks_for_all_segregative_sites(variants))
    print(time.time() - start)

if __name__ == "__main__":
    main()
