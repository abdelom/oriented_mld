import numpy as np
from detect_breakpoints import internal_incompatibility_2


def msprime_depth(tree_sequence):
    """

    :param tree_sequence: objet tree sequence de la suite tskit
    :return: liste des temps de coalescence, date des neouds du ts
    """
    return [node.time for node in tree_sequence.tables_dict["nodes"] if node.time > 0.0]


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
    list_incompatibilities = internal_incompatibilities(partition_dict, partition)
    for inf, sup in clade_in_block(list_incompatibilities, partition_dict[partition]):
        if inf is not None:
            if crown:
                delete = np.array(partition_dict[partition])
                delete = delete[delete > inf]
                delete = delete[delete < sup]
                delete = delete - inf
                variants_tmp = np.delete(variants[inf: sup + 1], delete)
                list_depth.append(depth_tmrca(variants_tmp, partition, len(variants[0].genotypes), len(partition) - 1))
            else:
                list_depth.append(depth_tmrca(variants[inf: sup], partition, len(variants[0].genotypes), len(partition)))
    return list_depth


def depths_clades(variants, crown, sample_size):
    """
    
    :param variants:  liste d'objet Variants du module msprime
    :param crown:  boolean, si true calcul du crown sinon du stems
    :param sample_size: nombre de qéuqence dans l'aligenement
    :return: ist des tmrca pour chaque bloque pour sous echantillon définie par la partition dans lequel la partition 
    corespondante est réalisée pour toute les partition rencontrée dans l'alignement
    """
    partition_dict = partition_position(variants, sample_size)
    list_depth = np.array([])
    for partition in partition_dict:
        list_depth = np.concatenate((list_depth, depths_clade(partition, partition_dict, variants, crown)))
    print(len(list_depth), len(list_depth[list_depth > 0]))
    return list_depth[list_depth > 0]


def depths_clades_size(variants, crown, sample_size, size):
    partition_dict = partition_position(variants, sample_size)
    list_depth_size = [np.array([])] * sample_size
    for partition in partition_dict:
        tmp = np.array(depths_clade(partition, partition_dict, variants, crown))
        list_depth_size[len(partition)] = np.concatenate((list_depth_size[len(partition)], tmp[tmp > 0]))
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

