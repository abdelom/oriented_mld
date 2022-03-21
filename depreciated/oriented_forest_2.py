import numpy as np
import multiprocessing as ml
from linckage import msprime_simulate_variants
import time

class Variant:
    def __init__(self, genotype):
        self.genotypes = genotype


def internal_incompatibility(genotype1, genotype2, oriented=True):
    A = {i for i in range(len(genotype1)) if genotype1[i]}
    B = {i for i in range(len(genotype2)) if genotype2[i]}
    return not (not A.intersection(B) or A.union(B) == B or A.union(B) == A)


def shortening(list_incompatible):
    list_incompatible = sorted(list_incompatible, key=lambda t: (t[0], t[1]))
    index = len(list_incompatible) - 1
    while index > 0:
        # print(list_incompatible)
        min1, max1 = list_incompatible[index - 1]
        min2, max2 = list_incompatible[index]
        if min1 <= min2 and max1 >= max2:
            list_incompatible[index - 1] = list_incompatible[index]
            del list_incompatible[index]
        elif min1 >= min2 and max1 <= max2:
            del list_incompatible[index]
        elif min1 < min2 and max1 <= max2 and max1 > min2:
            del list_incompatible[index]
            list_incompatible[index - 1] = (min2, max1)
        # elif min2 <= min1 and min2 >= max2 and max2 <= max1:
        #     del list_incompatible[index]
        #     list_incompatible[index - 1] = (min1, max2)
        index -= 1
    return list_incompatible


def detect_mld_block(i, variants, thresold=150):
    list_incompatible_sites = []
    j = i + 1
    while j < i + thresold and j < len(variants):
        if internal_incompatibility(variants[i].genotypes,
        variants[j].genotypes):
            list_incompatible_sites.append((i, j))
        j += 1
    return list_incompatible_sites


def mld_blocks(variants):
    start = time.time()
    pool = ml.Pool(4)
    data = pool.starmap(detect_mld_block, [(i, variants) for i in range(len(variants))])
    pool.close()
    print(time.time() - start)
    list_mld_blocks = [0] * len(variants)
    for start, end in shortening([item for sublist in data for item in sublist]):
        list_mld_blocks[start].append(end)
        list_mld_blocks[end].append(start)
    return list_mld_blocks


def rigth_breakpoint(variants, list_mld_blocks, index, segregated):
    i = index + 1
    while i < len(list_mld_blocks):

        i += 1

def left_breakpoint():
    pass


def mld_blocks_for_all_segregative_sites(list_mld_blocks, variants):
    sample_size = len(variants[0].genotypes)
    list_mld_for_segregatives = []
    for index, variant in enumerate(variants):
        tmp = sum(variant.genotypes)
        if tmp > 2 or tmp < sample_size:
            segregated = list(variant.genotypes == 1)
            nb_breakpoints_point = len(list_mld_blocks[index])
            if nb_breakpoints == 2:
                minimum = min(list_mld_blocks[index])
                maximum = max(list_mld_blocks[index])
                list_mld_for_segregatives.append((minimum, index), (index, maximum), segregated)
            elif nb_breakpoints == 1:
                if list_mld_blocks[index][0] < index:
                else:
            else:



def main():
    start = time.time()
    params = {"sample_size": 8, "Ne": 1, "ro": 8e-4, "mu": 8e-2,  "Tau": 1.0,
    "Kappa": 1.0 , "length": int(1e5)}
    variants = list(msprime_simulate_variants(params).variants())
    print(len(variants))
    print(len(mld_blocks(variants)))
    print(time.time() - start)


if __name__ == "__main__":
    main()
