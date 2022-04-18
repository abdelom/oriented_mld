import pandas as pd
import collections as cl
import matplotlib.pyplot as plt
from detect_breakpoints import *
from msprime_simulation import msprime_simulate_variants, test_tsinfer
from classification_breakpoints import class_brkpoints


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


def closest(obs_events, th_events, variants, length, chop=True):
    """
    retourne la moyennes des distances entre les événeents de recombinaisons et l'événement de recombinaison rélle le plus proche
    """
    positions, kind = zip(*obs_events)
    obs_events = np.array(positions)[np.array(kind) != "silent"]
    obs_events = [-length] + list(obs_events) + [2 * length]
    if chop:
        th_events = choping(th_events, variants, length)[1:-1]
    index = 0
    list_dist = []
    left = obs_events[0]
    for right in obs_events[1:]:
        while index < len(th_events) and left <= th_events[index] <= right:
            list_dist.append(min(th_events[index] - left, right - th_events[index]))
            index += 1
        left = right
    return np.mean(list_dist)


def th_events_ts_infer(obs_events, variants, params):
    """
    retourne le temps de calcul de tsinfer, le niombre de changements de topologies detecté et la distance au point de recombinaison
    discret ou incompatible le plus proche.
    """
    start = time.time()
    ts, edges, th_events, variants = test_tsinfer(variants, params["length"], simplify=True)
    end = time.time() - start
    dist_closest = closest(obs_events, th_events, variants, params["length"], chop=False)
    return [cl.Counter(list(zip(*class_brkpoints(edges, th_events, unroot=False)))[1])["incompatible"], \
            end, dist_closest]


def th_events_ag(obs_events, variants, params):
    """
    retourne le temps de calcul de tsinfer, le niombre de changements de topologies detecté et la distance au point de recombinaison
    discret ou incompatible le plus proche.
    """
    start = time.time()
    th_events = detect_events(variants, params["sample_size"])
    end = time.time() - start
    dist_closest = closest(obs_events, th_events, variants, params["length"])
    return len(th_events), end, dist_closest


def th_events_ek(obs_events, variants, params):
    """
    retourne le temps de calcul de la méhde naïve, le niombre de changements de topologies detecté et la distance au point de recombinaison
    discret ou incompatible le plus proche.
    """
    start = time.time()
    th_events = detect_internal_incompatibilities(variants, True, params["thresold"])
    end = time.time() - start
    dist_closest = closest(obs_events, th_events, variants, params["length"])
    return len(th_events), end, dist_closest


def data_simulation(params):
    closest_dist, result, list_time = [], [], []
    thresolds = [0]
    ts, edges, events, variants = msprime_simulate_variants(params)
    events = class_brkpoints(edges, events, unroot=True)
    events_count = cl.Counter(list(zip(*events))[1])
    total = events_count["incompatible"] + events_count["discret"]
    result.append(total)
    for func in [th_events_ts_infer, th_events_ag, th_events_ek]:
        if func == th_events_ek:
            thresolds = (50, 100)
        for thr in thresolds:
            params.update({"thresold": thr})
            nb_events, end, dist_closest = func(events, variants, params)
            closest_dist, result, list_time = closest_dist + [dist_closest], result + [nb_events], list_time + [end]
    return result, list_time, closest_dist


def method_comparaison(params, nb_repetition):
    i = 0
    list_times, list_results, list_closest = [], [], []
    list_result = []
    while i < nb_repetition:
        result, list_time, closest = data_simulation(params)
        print(result, list_time, closest)
        list_times.append(list_time.copy())
        list_result.append(result.copy())
        list_closest.append(closest.copy())
        i += 1
    pd.DataFrame(list_times, columns=["jerome", "abdel_guillaume", "elise_50", "elise_100"]).to_csv("time")
    pd.DataFrame(
        list_result, columns=["total", "jerome", "abdel_guillaume", "elise_50", "elise_100"]).to_csv("result")
    pd.DataFrame(
        list_closest, columns=["jerome", "abdel_guillaume", "elise_50", "elise_100"]).to_csv("closest")


def nb_partition(variants, individuals):
    dict_partitions = {}
    individuals = np.array(range(individuals))
    print(list(individuals))
    for variant in variants:
        partition = tuple(individuals[variant.genotypes == 1])
        if partition in dict_partitions:
            dict_partitions[partition] += 1
        else:
            dict_partitions[partition] = 0
    return dict_partitions


def scatter_plot(data):
    # data = pd.melt(results, id_vars=["total"],
    #                value_vars=["jerome", "abdel_guillaume", "elise_50", "elise_100"])
    fig, ax = plt.subplots()
    #colors = {'abdel_guillaume': 'red', 'elise_50': 'green', 'elise_100': 'blue', 'jerome': 'black'}
    plt.plot(data['total'], data['total'], color='blue', ls="-")
    ax.scatter(data['total'], data['jerome'], color="black")
    ax.scatter(data['total'], data['abdel_guillaume'], color="red")
    ax.scatter(data['total'], data['elise_50'], color="blue")
    ax.scatter(data['total'], data['elise_100'], color="green")
    plt.legend(['x = y', 'tsinfer', 'hierarchie', 'naive_50', 'naive_100'], title="Legend")
    fig.savefig("scatter.png", dpi=200)


def box_plot_time(data):
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.boxplot([data["jerome"], data["abdel_guillaume"], data["elise_50"],
                data["elise_100"]])  # , notch=True, patch_artist=True)
    plt.xticks([1, 2, 3, 4], ['tsinfer', 'hierarchie', 'naive_50', 'naive_100'])
    fig.savefig("box.png", dpi=200)


def covering(mld_list, length):
    cover = np.array([0] * length)
    inf = 0
    for tmp, sup in mld_list:
        cover[inf:sup] += 1
        inf = tmp
    cover[inf:length] += 1
    print(len(cover[cover == 3]) / len(cover) )
    fig, ax = plt.subplots(figsize=(15, 10))
    plt.plot(range(length), cover)
    fig.savefig("ddd.png", dpi=200)
    return cover


def main():
    params = {"sample_size": 10, "Ne": 1, "ro": 8e-6, "mu": 8e-4, "Tau": 1, "Kappa": 5, "length": int(1e7)}
    # method_comparaison(params, 300)
    # times, result = method_comparaison(params, 100)
    # times.to_csv("time")
    # result.to_csv("result")

    ts, edges, events, variants = msprime_simulate_variants(params)
    print(nb_partition(variants, 10))
    print(len(events))
    th_events = detect_events(variants, params["sample_size"])
    print(covering(th_events, int(1e7)))
    #obs_events = class_brkpoints(edges, events, unroot=True)
    #
    # results = pd.read_csv("result")
    # times = pd.read_csv("time")
    # closest = pd.read_csv("closest")
    # scatter_plot(results)
    # box_plot_time(times)
    # box_plot_time(closest)


if __name__ == "__main__":
    main()
