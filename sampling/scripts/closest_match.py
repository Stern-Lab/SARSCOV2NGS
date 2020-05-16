import argparse
from itertools import combinations
from Bio import Phylo
import numpy as np
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import os


def get_closest_cophetic_match(tree, tip_name, country):
    """
    find the closest match for each leave sequence belonging to country
    :param tree: a phylogenetic tree object
    :param tip: a tip in the tree to whom we want to find a match
    :param country: the country for whom the analysis is made
    :return: the label in the fasta file of the closest match
    """
    tip = [t for t in tree.get_terminals() if t.name == tip_name][0]
    node_path = tree.get_path(tip)
    parent_clade = node_path[-2]
    parent_clade_tips = [p for p in parent_clade.get_terminals() if country not in p.name.lower()] + [tip]

    # if homogeneous clade - keep going up until finding a mixed clade.
    i = -2
    limit = 10
    while len(parent_clade_tips) == len([t for t in parent_clade_tips if country in t.name.lower()]) or \
            len(parent_clade_tips) < limit:
        i -= 1
        parent_clade = node_path[i]
        parent_clade_tips = [p for p in parent_clade.get_terminals() if country not in p.name.lower()] + [tip]
    n = len(parent_clade_tips)

    # transform data to a dendrogram for cophenetic distance calculation
    # by - https://stackoverflow.com/questions/31033835/newick-tree-representation-to-scipy-cluster-hierarchy-linkage-matrix-format
    idx_dict = {}
    for i, t in enumerate(parent_clade_tips):
        idx_dict[t] = i
    label_dict = dict((str(y),x.name) for x,y in idx_dict.items())
    tip_idx = idx_dict[tip]
    dmat = np.zeros((n,n))
    for l1, l2 in combinations(parent_clade_tips, 2):
        d = tree.distance(l1, l2)
        dmat[idx_dict[l1], idx_dict[l2]] = dmat[idx_dict[l2], idx_dict[l1]] = d
    schlink = sch.linkage(squareform(dmat), method='average', metric='euclidean')

    # get distance between clusters. TODO check the benefit of this distance metric after tree evaluation.
    cophenet_distance = squareform(sch.cophenet(schlink))

    # extract the tip with the smallest cophenet distance
    # if all distances are the same - randomly pick one
    idx = np.argpartition(cophenet_distance[tip_idx], 1)[1]
    return label_dict[str(idx)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", type=str, required=True, help="path to a tree file")
    parser.add_argument("--country", type=str, default='Israel', help="country identifier")
    parser.add_argument("--name", type=str, required=True, help="tip identifier")
    parser.add_argument("--output", required=True, help="FASTA output file")
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    country = args.country.lower()

    match = get_closest_cophetic_match(tree, args.name, country)

    with open(args.output, "w") as output:
        output.write(match + '\n' + args.name)