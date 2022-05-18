import pandas as pd
import numpy as np
import collections as cl
import matplotlib.pyplot as plt
from detect_breakpoints import *
from msprime_simulation import msprime_simulate_variants, test_tsinfer
from tree_depth import *
from scipy.stats import kde


def events_counts(variants, obs_events, th_events):
    """
    retourn le nombre de point de recombinant contenue dans les intervalles définissant les jonctions entre les bloques
    """
    # nombre de points de recombinaisons des différents types dans chaque intervales délimitants les bloques mld
    dict_break_mld = {key: {"incompatible": 0, "silent": 0, "discret": 0, "hidden": 0} for key in th_events}
    for obs_event, event_type in obs_events:
        for inf, sup in th_events:
            if variants[inf].site.position >= obs_event or obs_event >= variants[sup].site.position:
                continue
            dict_break_mld[(inf, sup)][event_type] += 1
    return pd.DataFrame.from_dict(dict_break_mld, orient="index")  # cl.Counter(dict_break_mld.values())


def event_index(event, variants):
    index = 0
    while index < len(variants) and event > variants[index].site.position:
        index += 1
    return index + 1


def incompatibilities_in_block(list_nrc, variants, oriented):
    """
    retourne le nombre de bloques suposé non recombinant contenant des incompatibilités
    """
    count = 0
    for k, block in enumerate(list_nrc):
        print(k, len(list_nrc), count)
        index1, index2 = block
        for i, variant1 in enumerate(variants[index1:index2 + 1]):
            j = i + 1
            variant2 = variants[index1:index2 + 1][j]
            while j < len(variants[index1:index2 + 1]):
                variant2 = variants[index1:index2 + 1][j]
                j += 1
                if internal_incompatibility(variant1.genotypes, variant2.genotypes, oriented=oriented):
                    print(i, j, variant1.genotypes, variant2.genotypes, index1, index2)
                    count += 1
                    break
            if internal_incompatibility(variant1.genotypes, variant2.genotypes, oriented=oriented):
                break
    return count


def block_length(mld):
    """
    distributions des longueurs de blocs
    """
    length_distribution = [0] * 21
    lc = np.mean(mld)
    nb_blocks = len(mld)
    i = 0
    for length in mld:
        # print(lc, length, length / lc, (length / lc) // 0.1)
        index = (length / lc) // 0.2
        if index < 20:
            length_distribution[int(index)] += 1
        else:
            length_distribution[-1] += 1
        i += 1
    return np.array(length_distribution) / nb_blocks, lc


def plot_time_density(list_time, mu):
    fig, ax = plt.subplots(figsize=(15, 10))
    for elem in list_time:
        if list(elem):
            data = np.array(elem) / mu
            density = kde.gaussian_kde(data)
            x = np.linspace(0,  5, 2000)
            y = density(x)
            plt.plot(x, y)
    fig.savefig("time5.png", dpi=200)


def plot_time_scatter(list_time, dir):
    fig, ax = plt.subplots()
    colors = ["black", "red", "blue", "green"] * 4
    for index, elem in enumerate(list_time):
        #elem = block_length(elem)
        elem, lc = summary_stat(elem)
        print(elem)
        plt.plot(np.arange(0, len(elem) / 5, 0.2), elem, color=colors[index])
    fig.savefig(dir + ".png", dpi=200)


def main():
    params = {"sample_size": 6, "Ne": 1, "ro": 8e-6, "mu": 8e-4, "Tau": 1, "Kappa": 1, "length": int(1e8)}
    ts, edges, events, variants = msprime_simulate_variants(params)
    a = partition_dict_msprime(ts)
    b = partition_dict_estimate(variants, True, 6)
    start = time.time()
    msprime_estimate_scatters(estimate_v_msprime(b, a, 6), "a")
    print(time.time() - start)

    # times, result = method_comparaison(params, 100)
    # times.to_csv("time")
    # result.to_csv("result")
    list_time = []
    list_time_2 = []
    # for mu in range(-6, -3):
    #     params.update({"mu": 10 ** mu})
    #     tmp = []
    #     for i in range(10):
    #         ts, edges, events, variants = msprime_simulate_variants(params)
    #         tmp += list(depths_clades_size(variants, True, params["sample_size"], False))
    #     list_time.append(np.array(tmp))
    #     list_time_2 += msprime_depth(ts)
    # plot_time_scatter(list_time, "time_scatter_2")
    # for i in range(10):
    #     print(i)
    #     ts, edges, events, variants = msprime_simulate_variants(params)
    #     list_time += np.matrix([np.array(liste) for liste, _ in depths_clades_size(variants, True, params["sample_size"], True)])
    #     print(list_time)
    #     plot_time_scatter(list_time, "time_scatter_6", 0)
    # plot_time_density(list_time, params["mu"])
    list_time = []
    # for kappa in range(-1, 2):
    #     params.update({"Kappa": 10 ** kappa})
    #     ts, edges, events, variants = msprime_simulate_variants(params)
    #     list_time.append(np.array(depths_clades(variants, True, params["sample_size"])))
    #     list_time_2.append(msprime_depth(ts))
    # plot_time_scatter(list_time, "time_scatter_2")


if __name__ == "__main__":
    main()
