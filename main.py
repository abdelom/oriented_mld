import pandas as pd
import collections as cl
import tsinfer
import msprime
import tqdm
import itertools as it
from detect_breakpoints import *
from msprime_simulation import msprime_simulate_variants
from classification_breakpoints import class_brkpoints


def events_counts(variants: list, obs_events: list, th_events: list) -> list:
    dict_break_mld = {key: {"incompatible": 0, "silent": 0, "discret": 0, "hidden": 0} for key in th_events}
    for obs_event, event_type in obs_events:
        for inf, sup in th_events:
            if variants[inf].site.position >= obs_event or obs_event >= variants[sup].site.position:
                continue
            dict_break_mld[(inf, sup)][event_type] += 1
    return pd.DataFrame.from_dict(dict_break_mld, orient="index")  # cl.Counter(dict_break_mld.values())


def test_tsinfer(variants):
    with tsinfer.SampleData( sequence_length=int(1e6), num_flush_threads=2
    ) as sample_data:
        for var in variants:
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)
    ts = tsinfer.infer(sample_data)
    return ts.edges(), [brkpoint for brkpoint in ts.breakpoints()]


def incompatible_blocks_2(variants, oriented, list_mld):
    count = 0
    index1, index2 = 0, 1
    for br in list_mld[:-1]:
        variant1, variant2 = variants[index1], variants[index2]
        while variant1.site.position < br:
            while variant2.site.position < br:
                if internal_incompatibility(variant1.genotypes, variant2.genotypes, oriented=False):
                    count += 1
                    break
                index2 = index2 + 1
                variant2 = variants[index2]
            index1 += 1
            index2 = index1 + 1
            variant1 = variants[index1]
    return count


def main():
    params = {"sample_size": 10, "Ne": 1, "ro": 8e-5, "mu": 8e-4, "Tau": 1, "Kappa": 1, "length": int(1e6)}
    edges, breakpoints, variants = msprime_simulate_variants(params)
    th_events = detect_internal_incompatibilities(variants, oriented=True, thresold=50)
    obs_events = class_brkpoints(edges, breakpoints)
    data = events_counts(variants, obs_events, th_events)
    edges, br = test_tsinfer(variants)
    tmp = class_brkpoints(edges, br)
    print(incompatible_blocks_2(variants, True, br))
    # nrc_blocks = [(0, 0)] + th_events + [(len(variants) - 1, 0)]
    # nrc_blocks = [(nrc_blocks[i - 1][1], nrc_blocks[i][0]) for i in range(1, len(nrc_blocks))]
    # data_2 = events_counts(variants, obs_events, nrc_blocks)
    # print(data_2)
    # for i in range(len(data["breakpoint"])):
    # print(tuple(data.iloc[i]))


if __name__ == "__main__":
    main()
