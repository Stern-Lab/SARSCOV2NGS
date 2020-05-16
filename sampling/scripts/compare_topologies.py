import argparse
from itertools import combinations
from scipy.spatial import distance
import pandas as pd
import numpy as np
import os
from Bio import SeqIO, Phylo

def prune_world_tree(tree, sampled_tree):
    """
    prune the world tree to include only tips from a sampled tree
    :param tree: the complete world tree
    :param sampled_tree: the down-sampled tree
    :return: the world tree after tips prunning
    """
    sampled_tips = [t.name for t in sampled_tree.get_terminals()]
    prune_tips = [t for t in tree.get_terminals() if t.name not in sampled_tips]
    for t in prune_tips:
        tree.prune(t)
    return tree


def map_tree_to_vector_idx(tree):
    """
    map all pairs of tips of a tree to their corresponding index in the tree vector
    :param tree: a ML tree
    :return: a dictionary where (t1,t2) are sorted keys and the value is an index.
    """
    mapper = {}
    tips = [t.name for t in tree.get_terminals()]
    i = 0
    for t1, t2 in combinations(tips, 2):
        name = tuple(sorted([t1, t2]))
        mapper[name] = i
        i += 1
    return mapper


def tree_to_vector(tree, tip_mapping):
    """
    implement the Kendall-Colijin metric
    :param tree: a ML tree
    :param tip_mapping: mapping of each tip to a vector index
    :return: a vector representation of the tree
    """
    tips = tree.get_terminals()
    # init vector
    vec = [0] * len(tip_mapping)
    for t1, t2 in combinations(tips, 2):
        mrca = tree.common_ancestor(t1,t2)
        root_2_mrca_distance = len(tree.get_path(mrca))
        name = tuple(sorted([t1.name, t2.name]))
        vec[tip_mapping[name]] = root_2_mrca_distance

    # add 1's at the end of the vector for each tip
    vec = vec + [1] * len(tips)
    return vec

def get_trees_distance(tree1, tree2, tree_mapping):
    """
    given two trees tree1, tree2 return their distance according to Kendall- Colijin metric
    :param tree1: first tree
    :param tree2: second tree
    :param tree_mapping: mapping of tip pairs to their vector index
    :return: topology tree distance
    """
    vec1 = tree_to_vector(tree1, tree_mapping)
    vec2 = tree_to_vector(tree2, tree_mapping)
    dist = distance.euclidean(vec1, vec2)
    return dist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--world_tree", type=str, required=True, help="path to a tree file")
    parser.add_argument("--sampled_tree", type=str, required=True, help="path to a tree file")
    parser.add_argument("--output", required=True, help="FASTA output file")
    args = parser.parse_args()

    tree1 = Phylo.read(args.world_tree, 'newick')
    tree2 = Phylo.read(args.sampled_tree, 'newick')

    tree1.root_with_outgroup({'name': 'Wuhan-Hu-1/2019'})
    tree2.root_with_outgroup({'name': 'Wuhan-Hu-1/2019'})

    tree1 = prune_world_tree(tree1, tree2)

    mapper = map_tree_to_vector_idx(tree2)
    dist = get_trees_distance(tree1, tree2, mapper)
    np.savetxt(args.output, dist)




