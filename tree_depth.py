import numpy as np
from detect_breakpoints import internal_incompatibility_2, detect_events
import matplotlib.pyplot as plt
from msprime_simulation import msprime_simulate_variants, test_tsinfer
from sklearn.linear_model import LinearRegression


def rsquared(data, dict_reg, dict_pred=None):
    for key in data:
        reg = LinearRegression(fit_intercept=False).fit(data[key][0].reshape(len(data[key][0]), 1),
                                                        data[key][1].reshape(len(data[key][0]), 1))
        dict_reg[key].append(reg.score(data[key][0].reshape(len(data[key][0]), 1), data[key][1].reshape(len(data[key][0]), 1)))
        if dict_pred:
            dict_pred[key] = reg.predict(data[key][0].reshape(len(data[key][0]), 1))
    return dict_reg, dict_pred


def msprime_depth(tree_sequence):
    """

    :param tree_sequence: objet tree sequence de la suite tskit
    :return: liste des temps de coalescence, date des neouds du ts
    """
    return [node.time for node in tree_sequence.tables_dict["nodes"] if node.time > 0.0]

########################################################################################################################
##########################################estimateur##################################################################
########################################################################################################################

def depth_tmrca(variants, partition, sample_size, k):
    """

    :param variants: list d'objets Variants de la suite tskit, liste des sites polymorphes
    :param partition: sous échantillon des individues de l'aligenement
    :param sample_size: le nombre d'individues total
    :param k: la taille du sous échantillon, pour le stem, ou la taille du sous échantillon - 1, pour l crown
    :return: depth, float, la profondeur, au taux de mutation près, de du tmrca au sous échantillon
    """
    depth = 0
    length = variants[-1].site.position - variants[0].site.position
    individuals = np.array(range(sample_size))
    partition = set(partition)
    for variant in variants[1:-1]:
        partition_tmp = set(individuals[variant.genotypes == 1])
        if partition_tmp.union(partition) == partition:
            depth += len(partition_tmp)
    return depth / (length * k)


def crown_stem(crown, inf, sup, variants, partition, partition_pos):
    if crown:
        delete = np.array(partition_pos)
        delete = delete[delete > inf]
        delete = delete[delete < sup]
        delete = delete - inf
        variants_tmp = np.delete(variants[inf: sup + 1], delete)
        return depth_tmrca(variants_tmp, partition, len(variants[0].genotypes), len(partition))
    else:
        return depth_tmrca(variants[inf: sup], partition, len(variants[0].genotypes), len(partition))

########################################################################################################################
##########################################hierarchical bloc#############################################################
########################################################################################################################


########################################################################################################################
##########################################maximal bloc##################################################################
########################################################################################################################

def partition_position(variants, sample_size):
    """

    :param variants: list d'objets Variants de la suite tskit, liste des sites polymorphes
    :param sample_size: le nombre d'individues total
    :return:partition_dict: dictionnaire, les clé sont les partitions engendrés par les snps et les valeurs du dictinnaire
    les listes des snp générant les partition
    """
    individuals = np.array(range(sample_size))
    partition_dict = {}
    for index, variant in enumerate(variants):
        partition = tuple(individuals[variant.genotypes == 1])
        if len(partition) < 2:
            continue
        if partition in partition_dict:
            partition_dict[partition].append(index)
        else:
            partition_dict[partition] = [index]
    return partition_dict


def internal_incompatibilities(partition_dict, partition):
    """

    :param partition_dict: dictionnaire, les clé sont les partitions engendrés par les snps et les valeurs du dictinnaire
    :param partition: une partition corespondantz à un snp rencontré au moins une fois dans l'alignement
    :return: la liste des positions de tout les snp incompatibles avec la partition courante
    """
    list_incompatibilities = []
    for partition_2 in partition_dict:
        if internal_incompatibility_2(set(partition), set(partition_2)): # si les
            list_incompatibilities += partition_dict[partition_2]
    return sorted(list_incompatibilities)


def clade_in_block(list_incompatibilities, partition_pos):
    """

    :param list_incompatibilities:liste des sites incompatible avec la ti courante.
    :param partition_pos: liste des position de la partition courante
    :return:list de tuples délimitant des bloques dans lequel la partition est réalisée
    """
    list_blocks = []
    # print(list_incompatibilities)
    inf = list_incompatibilities[0]
    i = 0
    for sup in list_incompatibilities[1:]:
        while i < len(partition_pos) and partition_pos[i] < inf:
            i += 1
        if i == len(partition_pos):
            break
        if inf < partition_pos[i] < sup:
            list_blocks.append((inf, sup))
        inf = sup
    return list_blocks


def depths_clade(partition, partition_dict, variants, crown):
    """

    :param partition: partition courante
    :param partition_dict: dictionnaire, les clé sont les partitions engendrés par les snps et les valeurs du dictinnaire
    :param variants:  liste d'objet Variants du module msprime
    :param crown: boolean, si true calcul du crown sinon du stems
    :return: list des tmrca pour chaque bloque pour sous echantillon définie par la partition dans klequel la partition correspondante est réalisé
    """
    list_depth = []
    list_born = []
    list_incompatibilities = [elem for elem in internal_incompatibilities(partition_dict, partition) if elem]
    for inf, sup in clade_in_block(list_incompatibilities, partition_dict[partition]):
        if inf is not None:
            if crown:
                delete = np.array(partition_dict[partition])
                delete = delete[delete > inf]
                delete = delete[delete < sup]
                delete = delete - inf
                variants_tmp = np.delete(variants[inf: sup + 1], delete)
                list_depth.append(depth_tmrca(variants_tmp, partition, len(variants[0].genotypes), len(partition)))
            else:
                list_depth.append(depth_tmrca(variants[inf: sup], partition, len(variants[0].genotypes), len(partition)))
            list_born.append((inf, sup))
    list_depth = np.array(list_depth)
    return np.array(list_born)[list_depth > 0], list_depth[list_depth > 0]


def depths_clades_size(variants, crown, sample_size, size):
    partition_dict = partition_position(variants, sample_size)
    list_depth_size = [np.array([])] * sample_size
    for partition in partition_dict:
        tmp = np.array(depths_clade(partition, partition_dict, variants, crown)[1])
        list_depth_size[len(partition)] = np.concatenate((list_depth_size[len(partition)], tmp))
    if size:
        list_stat = []
        for list_depth in list_depth_size:
            if not list(list_depth):
                continue
            list_stat.append(summary_stat(list_depth))
        return list_stat
    return np.concatenate(list_depth_size)


def summary_stat(list_h):
    """
    distributions des longueurs de blocs
    """
    length_distribution = [0] * 21
    lc = np.mean(list_h)
    nb_blocks = len(list_h)
    i = 0
    for length in list_h:
        # print(lc, length, length / lc, (length / lc) // 0.1)
        index = (length / lc) // 0.2
        if index < 20:
            length_distribution[int(index)] += 1
        else:
            length_distribution[-1] += 1
        i += 1
    return np.array(length_distribution) / nb_blocks, lc


########################################################################################################################
##########################################msprime vs estimate###########################################################
########################################################################################################################

def clade_time_msprime(tree, root):
    if tree.is_leaf(root):
        return [([root], {"time": 0})]
    clade_left = clade_time_msprime(tree, tree.left_child(root))
    clade_right = clade_time_msprime(tree, tree.right_child(root))
    new_clade = [(clade_left[0][0] + clade_right[0][0], {"time": tree.time(root), "node": root})]
    return new_clade + clade_right + clade_left


def partition_dict_msprime(ts):
    hierarchie = {}
    for tree in ts.trees():
        left, right = tree.interval.left, tree.interval.right
        clade = {tuple(sorted(key)): value for key, value in clade_time_msprime(tree, tree.root) if value["time"] > 0}
        for partition in clade:
            if partition in hierarchie:
                if left == hierarchie[partition][-1]["right"] and hierarchie[partition][-1]["node"] == clade[partition]["node"]:
                    hierarchie[partition][-1].update({"right": right})
                else:
                    clade[partition].update({"left": left, "right": right})
                    hierarchie[partition].append(clade[partition])
            else:
                clade[partition].update({"left": left, "right": right})
                hierarchie[partition] = [clade[partition]]
    return hierarchie


def partition_dict_tsinfer():
    pass


def partition_dict_estimate(variants, crown, sample_size):
    partition_dict = partition_position(variants, sample_size)
    dict_depth_size = {}
    for partition in partition_dict:
        borns, depths = depths_clade(partition, partition_dict, variants, crown)
        for depth, born in zip(depths, borns):
            if partition in dict_depth_size:
                dict_depth_size[partition].append({"time": depth, "right": variants[born[1]].site.position, "left":  variants[born[0]].site.position})
            else:
                dict_depth_size[partition] = [{"time": depth, "right":  variants[born[1]].site.position, "left":  variants[born[0]].site.position}]
    return dict_depth_size


def depht_hierarchie(crown, dict_depth_size, variants, sample_size):
    partition_dict = partition_position(variants, sample_size)
    for partition in dict_depth_size:
        for dictionary in dict_depth_size[partition]:
            inf, sup = dictionary["left"], dictionary["right"]
            dictionary.update({"time": crown_stem(crown, inf, sup, variants, partition, partition_dict[partition]),
                               "left": variants[inf].site.position,
                               "right": variants[sup].site.position})
    return dict_depth_size


def partition_dict_hierarchie(variants, crown, sample_size):
    list_blocks = [(0, 0)] + detect_events(variants, sample_size)
    list_blocks.append((len(list_blocks) - 1, 0))
    dict_depth_size = {}
    individuals = np.array(range(sample_size))
    for inf, sup in [(list_blocks[i][0], list_blocks[i + 1][1]) for i in range(len(list_blocks) - 1)]:
        for variant in variants[inf:sup]:
            partition = tuple(individuals[variant.genotypes == 1])
            if partition in dict_depth_size:
                if inf < dict_depth_size[partition][-1]["right"]:
                    dict_depth_size[partition][-1].update({"right": sup})
                else:
                    dict_depth_size[partition].append({"left": inf, "right": sup})
            elif len(partition) > 1:
                dict_depth_size[partition] = [{"left": inf, "right": sup}]
    return depht_hierarchie(crown, dict_depth_size, variants, sample_size)


def partition_estimate_msprime(list_estimate, list_msprime, mu):
    estimator, ms = [], []
    index = 0
    for estimate in list_estimate:
        while index < len(list_msprime):
            mspr = list_msprime[index]
            if mspr["right"] < estimate["right"] and mspr["left"] > estimate["left"]:
                estimator = np.append(estimator, estimate["time"])
                ms = np.append(ms, mspr["time"] * mu)
            if mspr["right"] > estimate["right"]:
                break
            index += 1
    return estimator, ms


def estimate_v_msprime(estimate_msprime, dict_estimate, dict_msprime, mu):
    for partition in dict_estimate:
        estimate, ms = partition_estimate_msprime(dict_estimate[partition], dict_msprime[partition], mu)
        estimate_msprime[len(partition)][0] = np.concatenate([estimate_msprime[len(partition)][0], ms])
        estimate_msprime[len(partition)][1] = np.concatenate([estimate_msprime[len(partition)][1], estimate])
    return estimate_msprime


def replicate_ems(params, nb):
    dict_rand = {}
    for mu_ro in [1, 10, 100]:
        estimate_msprime = dict([(i + 1, [np.array([]), np.array([])]) for i in range(1, params["sample_size"] - 1)])
        params.update({"mu": mu_ro * 8e-6})
        for i in range(nb):
            ts, _, _, variants = msprime_simulate_variants(params)
            a = partition_dict_msprime(ts)
            b = partition_dict_hierarchie(variants, True, 6)
            estimate_msprime = estimate_v_msprime(estimate_msprime, b, a, params["mu"])
            dict_rand[mu_ro] = estimate_msprime
        msprime_estimate_scatters(estimate_msprime, mu_ro, "fig_dir")
        cumulative_times(estimate_msprime, mu_ro, "fig_dir")
    residues_box(dict_rand, mu_ro, "fig_dir",  params["sample_size"])


def comp_method(dict_method, params, nb):
    ts, _, _, variants = msprime_simulate_variants(params)
    dict_method_r = {}
    for key in dict_method:
        dict_rsquared = {i: [] for i in range(2, params["sample_size"])}
        for i in range(nb):
            estimate_msprime = dict([(i + 1, [np.array([]), np.array([])]) for i in range(1, params["sample_size"] - 1)])
            ts, _, _, variants = msprime_simulate_variants(params)
            estimate_msprime = estimate_v_msprime(estimate_msprime, dict_method[key](variants, True, 6),
                                                  partition_dict_msprime(ts), params["mu"])
            dict_rsquared, _ = rsquared(estimate_msprime, dict_rsquared)
        dict_method_r[key] = dict_rsquared
    print(dict_method_r)
    return dict_method_r


########################################################################################################################
##########################################plot###########################################################
########################################################################################################################


def msprime_estimate_scatters(dict_times, mu_ro, dir_fig):
    fig, axs = plt.subplots(2, 2)
    for key in dict_times:
        data = dict_times[key]
        r2, prediction = rsquared(data)
        axs[(key - 2) // 2, (key - 2) % 2].plot(data[0], data[0], color='blue', ls="-")
        axs[(key - 2) // 2, (key - 2) % 2].plot(data[0], prediction, color='red', ls="-")
        axs[(key - 2) // 2, (key - 2) % 2].scatter(data[0], data[1], alpha=0.1)
        axs[(key - 2) // 2, (key - 2) % 2].set_title(f'k={key}, r2={round(r2, 2)}')
    for ax in axs.flat:
        ax.set(xlabel='msprime time', ylabel='estimate times')
    for ax in axs.flat:
        ax.label_outer()
    fig.savefig(f"{dir_fig}/scatter_{mu_ro}.png", dpi=200)


def rsquared_box(dict_squared, mu_ro, dir_fig, sample_size):
    fig, axs = plt.subplots(2, 2)
    labels = dict_squared.keys()
    for effective in range(2, sample_size):
        data = []
        for key in dict_squared:
            data.append(dict_squared[key][effective])
        axs[(effective - 2) // 2, (effective - 2) % 2].boxplot(data)
        #axs[(effective - 2) // 2, (effective - 2) % 2].set_xticks([1, 2, 3], labels)
        for ax in axs.flat:
            ax.set(xlabel='msprime time', ylabel='estimate times')
        for ax in axs.flat:
            ax.label_outer()
    fig.savefig(f"{dir_fig}/rsquared_{mu_ro}.png", dpi=200)


def cumulative_times(dict_times, mu_ro, dir_fig):
    fig, axs = plt.subplots()
    key = 2
    #for key in dict_times:
    data = sorted(dict_times[key][1])
    plt.plot(data, np.arange(0, 1,  1 / len(data)), color='blue', ls="-")
    fig.savefig(f"{dir_fig}/cumulative_{mu_ro}.png", dpi=200)


def residues_box(dict_times, mu_ro, dir_fig, sample_size):
    fig, axs = plt.subplots(2, 2)
    for effective in range(2, sample_size):
        data = []
        for key in dict_times:
            data.append((dict_times[key][effective][1] - dict_times[key][effective][0]) / dict_times[key][effective][0])
        axs[(effective - 2) // 2, (effective - 2) % 2].boxplot(data)
        #axs[(effective - 2) // 2, (effective - 2) % 2].set_xticks([1, 2, 3])
        for ax in axs.flat:
            ax.set(xlabel='msprime time', ylabel='estimate times')
        for ax in axs.flat:
            ax.label_outer()
    fig.savefig(f"{dir_fig}/sbox_residue_{mu_ro}.png", dpi=200)


def box_plot_time(data, mu_ro, out_dir):
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.boxplot([data["jerome"], data["abdel_guillaume"], data["elise_50"],
                data["elise_100"]])  # , notch=True, patch_artist=True)
    plt.xticks([1, 2, 3, 4], ['tsinfer', 'hierarchie', 'naive_50', 'naive_100'])
    fig.savefig(f"{out_dir}/box_{mu_ro}.png", dpi=600)


def main():
    params = {"sample_size": 6, "Ne": 1, "ro": 8e-6, "mu": 8e-5, "Tau": 1, "Kappa": 1, "length": int(1e8)}
    out_csv, out_fig = "csv_dir", "fig_dir"
    #replicate_ems(params, 1)
    rsquared_box(comp_method({"hierarchies": partition_dict_hierarchie, "estimate": partition_dict_estimate}, params, 30), 1, out_fig, 6)


if __name__ == "__main__":
    main()

