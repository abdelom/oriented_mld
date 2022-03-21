# import multiprocessing as ml
# import time
# import numpy as np
import pandas as pd
import collections as cl
import itertools as it
from detect_breakpoints import *
from msprime_simulation import msprime_simulate_variants
from classification_breakpoints import class_brkpoints


def recombination_events(params):
    edges, breakpoints, variants = msprime_simulate_variants(params)
    breakpoint_type = class_brkpoints(edges, breakpoints)
    th_events = detect_internal_incompatibilities(variants, oriented=False)
    detected_events, index = [-1] * (len(breakpoints) - 2), 0
    breakpoints = breakpoints[1:-1]
    inf, sup = th_events[0]
    for index_2, event in enumerate(breakpoints):
        while event > variants[sup].site.position and index < len(th_events) - 1:
            index += 1
            inf, sup = th_events[index]
        if variants[inf].site.position < event < variants[sup].site.position:
            detected_events[index_2] = (inf, sup)
    return pd.DataFrame({"breakpoint": breakpoints,
                         "type": breakpoint_type,
                         "detection": detected_events})


def accuracy(variants, breakpoints, th=20, oriented=True):
    mld = detect_internal_incompatibilities(variants, oriented, th)
    # mld = incompatible_blocks_2(variants, oriented, mld)
    nb_recombinaison = len(breakpoints)
    count = 0
    for brkpoint in breakpoints:
        for mini, maxi in mld:
            if variants[mini].site.position < brkpoint < variants[maxi].site.position:
                count += 1
                break
    return count / nb_recombinaison


def breakpoint_by_incompatility(variants, breakpoints, mld):
    dict_break_mld = {key: 0 for key in mld}
    for brkpoint in breakpoints:
        for min, max in mld:
            if variants[min].site.position < brkpoint < variants[max].site.position:
                dict_break_mld[(min, max)] += 1
                # mld.remove((min, max))
                # break
    return cl.Counter(dict_break_mld.values())


def detection_frequency(params, oriented=True):
    list_ro = [params["mu"] * 10 ** power for power in range(-3, 2, 1)]
    detected_events = {ro: 0 for ro in list_ro}
    for ro in list_ro:
        params.update({'ro': ro})
        _, breakpoints, variants = msprime_simulate_variants(params)
        mld = detect_internal_incompatibilities(variants, oriented=oriented, thresold=20)
        breackpoints_mld = breakpoint_by_incompatility(variants, breakpoints, mld=mld)
        breakpoints = len(breakpoints)
        detected_events[ro] += sum(np.array(list(breackpoints_mld.keys())) *
            np.array(list(breackpoints_mld.values()))) / breakpoints
    return detected_events


def block_length(mld):
    mld_length = []
    for i in range(1, len(mld)):
        if mld[i] - mld[i - 1] < 0:
            print(mld[i], mld[i - 1])
        mld_length.append(mld[i] - mld[i - 1])
    # mld_length = [mld[i] - mld[i - 1] for i in range(1, len(mld))]
    # print(mld[-4:])
    length_distribution = [0] * 26
    lc = np.mean(mld_length)
    nb_blocks = len(mld_length)
    i = 0
    for length in mld_length:
        # print(lc, length, length / lc, (length / lc) // 0.1)
        index = (length / lc) // 0.2
        if index < 25:
            length_distribution[int(index)] += 1
        else:
            length_distribution[-1] += 1
        i += 1
    return np.array(length_distribution) / nb_blocks


def incompatible_blocks_2(variants, oriented, list_mld):
    liste = []
    for i in range(len(list_mld) - 1):
        b_inf = list_mld[i][1]
        b_sup = list_mld[i + 1][0]
        for variant1, variant2 in it.combinations(variants[b_inf: b_sup], 2):
            if internal_incompatibility(variant1.genotypes, variant2.genotypes, oriented=oriented):
                tmp = detect_internal_incompatibilities(variants[b_inf: b_sup], thresold=len(variants[b_inf: b_sup]))
                tmp2 = [(elem[0] + b_inf, elem[1] + b_inf) for elem in tmp]
                liste += tmp2
                break
    return sorted(list_mld + liste, key=lambda t: (t[0], t[1]))


# def filter_incompatibiliies(incompatibilities, segregated, index, variants, oriented):
#     incompatibilities_2 = []
#     for incompatibility in incompatibilities:
#         index_1, index_2 = incompatibility
#         if internal_incompatibility(vrariants[index_1].genotypes[segregated], \
#         variants[index_2].genotypes[segregated], oriented=oriented) or index in incompatibility:
#             incompatibilities_2.append(incompatibility)
#     return shortening(incompatibilities_2)
#
# def find_block(list_mld, segregated, variants)
#     for index, variant in enumerate(variants):
#         if sum(variant.genotypes) >= 3:
#             segregated = (variants.genotypes == 1)
#             tmp = filter_incompatibiliies(list_mld, segregated, index, variants, True)
#             for elem in tmp:
#
#     return None

# def faux_positifs(params):
#     ts, breakpoints, variants = msprime_simulate_variants(params)
#     tmp, short = detect_internal_incompatibilities(variants, oriented=True, thresold=20)
#     #t = breakpoint_by_incompatility(variants, breakpoints, incompatible_sites, th=20)
#     v = breakpoint_by_incompatility(variants, breakpoints, short)
#     liste = []
#     for key, elem in v.items():
#         # print(key, elem)
#         if elem ==0:
#             liste.append(key)
#     return liste, variants


def main():
    params = {"sample_size": 10, "Ne": 1, "ro": 8e-5, "mu": 8e-4, "Tau": 1, "Kappa": 1, "length": int(1e6)}
    # ts, breakpoints, variants = msprime_simulate_variants(params)
    # print(class_brkpoints(ts.edges(), breakpoints))
    data = recombination_events(params)
    for i in range(len(data["breakpoint"])):
        print(tuple(data.iloc[i]))

if __name__ == "__main__":
    main()
