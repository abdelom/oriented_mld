import numpy as np
from detect_breakpoints import internal_incompatibility_2


def msprime_depth(tree_sequence):
    """

    :param tree_sequence: objet tree sequence de la suite tskit
    :return: liste des temps de coalescence, date des neouds du ts
    """
    """
    :param tree_sequence: 
    :return: 
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
        partition_tmp = set(individuals[variant.genotypes])
        if partition_tmp.union(partition) == partition:
            depth += len(partition_tmp) / (length * k)
    return depth


def partition_position(variants, sample_size):
    """

    :param variants: list d'objets Variants de la suite tskit, liste des sites polymorphes
    :param sample_size: le nombre d'individues total
    :return:partition_dict:
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

    :param partition_dict:
    :param partition:
    :return:
    """
    list_incompatibilities = []
    for partition_2 in partition_dict:
        if internal_incompatibility_2(set(partition), set(partition_2)):
            list_incompatibilities += partition_dict[partition_2]
    return sorted(list_incompatibilities)


def clade_in_block(list_incompatibilities, partition_pos):
    """

    :param list_incompatibilities:
    :param partition_pos:
    :return:
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

    :param partition:
    :param partition_dict:
    :param variants:
    :param crown:
    :return:
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
                list_depth.append(depth_tmrca(variants_tmp, partition, len(partition), len(partition) - 1))
            else:
                list_depth.append(depth_tmrca(variants[inf, sup], partition, len(partition), len(partition)))
    return list_depth


def depths_clades(variants, crown, sample_size):
    """
    
    :param variants:
    :param crown:
    :param sample_size:
    :return:
    """
    partition_dict = partition_position(variants, sample_size)
    list_depth = np.array([])
    for partition in partition_dict:
        list_depth = np.concatenate((list_depth, depths_clade(partition, partition_dict, variants, crown)))
    return list_depth[list_depth > 0]
