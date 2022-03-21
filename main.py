import pandas as pd
import collections as cl
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


def main():
    params = {"sample_size": 10, "Ne": 1, "ro": 8e-5, "mu": 8e-4, "Tau": 1, "Kappa": 1, "length": int(1e6)}
    edges, breakpoints, variants = msprime_simulate_variants(params)
    th_events = detect_internal_incompatibilities(variants, oriented=True, thresold=50)
    obs_events = class_brkpoints(edges, breakpoints)
    data = events_counts(variants, obs_events, th_events)
    print(data)
    nrc_blocks = [(0, 0)] + th_events + [(len(variants) - 1, 0)]
    nrc_blocks = [(nrc_blocks[i - 1][1], nrc_blocks[i][0]) for i in range(1, len(nrc_blocks))]
    data_2 = events_counts(variants, obs_events, nrc_blocks)
    print(data_2)
    # for i in range(len(data["breakpoint"])):
    # print(tuple(data.iloc[i]))


if __name__ == "__main__":
    main()
