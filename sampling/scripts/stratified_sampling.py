import argparse
from itertools import combinations
from Bio import SeqIO, Phylo
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
from tqdm import tqdm
import hashlib
import os


def get_closest_cophetic_match(tree, tip, country):
    """
    find the closest match for each leave sequence belonging to country
    :param tree: a phylogenetic tree object
    :param tip: a tip in the tree to whom we want to find a match
    :param country: the country for whom the analysis is made
    :return: the label in the fasta file of the closest match
    """
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



# TODO manual version (based on densities) - need to make it more generic
def bin_metadata_by_week(metadata):
    """
    bin metadata into time intervals for sampling
    :param metadata: metadata table
    :return: the same table with bin labels.
    """
    weeks = [g for n, g in metadata.groupby(pd.Grouper(key='datetime', freq='W'))]
    # combine the first 10 weeks of the epidemic into one
    bins = [pd.concat(weeks[:10]), weeks[10], weeks[11], weeks[12], pd.concat(weeks[13:])]
    for i, bin in enumerate(bins):
        bin['bin'] = i
    return pd.concat(bins)


def time_sampling(sequences, metadata, country, num_sequences):
    """
    stratified sample by time of the sequences, according to their dates in the metadata file
    """
    strains = [rec.id for rec in SeqIO.parse(sequences, 'fasta')]   # filter the data by the latest sequences filter
    metadata = metadata[(metadata['strain'].isin(strains)) | (metadata['country'] == country)]
    intervals = 5
    num_samples = num_sequences // intervals
    metadata['datetime'] = metadata['date'].apply(lambda x: pd.to_datetime(x, format="%Y-%m-%d"))

    # bin by weeks
    binned = bin_metadata_by_week(metadata)
    sampled = []
    for i in tqdm(binned['bin'].unique()):
        chosen = list(np.random.choice(binned[binned['bin'] == i]['strain'], num_samples))
        sampled.extend(chosen)
    # add all the sequences from the country to the sampled batch
    sampled = sampled + [s for s in strains if country in s]
    return sampled


def closest_match_sampling(tree, country):
    """
    get the closest match for each tip in the tree belonging to country
    """
    tips = [t for t in tree.get_terminals() if country in t.name.lower()]
    sampled = []
    for t in tqdm(tips):
        match = get_closest_cophetic_match(tree, t, country)
        sampled.append(match)
    return sampled


def sequence_2_md5(records, name):
    curr_rec = [n for n in records if n.name == name]
    if len(curr_rec) == 0:
        return 'invalid'
    else:
        return hashlib.md5(str(curr_rec[0].seq).encode('utf-8')).hexdigest()


def find_redundency(metadata, sequences):
    records = [r for r in SeqIO.parse(sequences, 'fasta')]
    isr = metadata[metadata['country'] == 'Israel']
    isr['seq_id'] = isr['strain'].apply(lambda x: sequence_2_md5(records, x))
    isr = isr[isr['seq_id'] != 'invalid']
    grouped = isr.groupby(['seq_id', 'originating_lab']).size().reset_index(name='size')

    seqs = isr[isr['seq_id'].isin(grouped[grouped['size'] > 1]['seq_id'].values)]\
    [['strain', 'originating_lab', 'date']]


    return seqs

# def closest_match_cluster():
#     tree = '/sternadi/nobackup/volume1/covid/ncov_ISR_most_updated/results/tree_raw.nwk'
#     tree_obj = Phylo.read(tree, 'newick')
#     tip_names = [t.name for t in tree_obj.get_terminals() if 'Israel' in t.name]
#     out_base = '/sternadi/home/volume3/COVID19/data/sampling/tip_matches/'
#     for t in tip_names:
#         script_runner(
#             'python /sternadi/home/volume3/COVID19/code/sampling/closest_match.py --tree {} --name {} --output {}'.format(
#                 tree, t, os.path.join(out_base, t.replace('/', '_') + '_match.txt')), alias='closest_match_covid',
#             load_python=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences after filtration")
    parser.add_argument("--exogenous", type = int, default=50, help="number of exogenous sequences to sample")
    parser.add_argument("--metadata", type=str, required=True, help="path to a metadata file")
    parser.add_argument("--tree", type=str, required=True, help="path to a tree file")
    parser.add_argument("--country", type=str, default='Israel', help="country identifier")
    parser.add_argument("--cluster", type=int, default=0, help="flag for tip values previously calculated")
    parser.add_argument("--output", required=True, help="FASTA output file")
    args = parser.parse_args()

    metadata = pd.read_table(args.metadata)
    tree = Phylo.read(args.tree, "newick")
    # make sure the tree is rooted
    tree.root_with_outgroup({'name': 'Wuhan-Hu-1/2019'})
    country = args.country.lower()
    metadata['country'] = metadata['country'].apply(lambda x: x.lower())

    sampled_strains = time_sampling(args.sequences, metadata, country, args.exogenous)
    # add to the list the strains from closest match
    if args.cluster == 0:
        with open(r'/Volumes/STERNADILABHOME$/volume3/COVID19/data/sampling/include.txt','r') as f:
            sampled_matches = [s.strip() for s in f.readlines()]
        with open(r'/Volumes/STERNADILABHOME$/volume3/COVID19/data/sampling/exclude.txt','r') as f:
            exclude = [s.strip() for s in f.readlines()]
        additional_tips = [t.name for t in tree.get_terminals() if country in t.name.lower()]
        sampled_matches = sampled_matches + [t for t in additional_tips if t not in exclude]
    else:
        sampled_matches = closest_match_sampling(tree, country)

    sampled = list(set(sampled_strains + sampled_matches + ['Wuhan-Hu-1/2019']))    # add the root

    records = [rec for rec in SeqIO.parse(args.sequences, 'fasta') if rec.id in sampled]
    cnt = SeqIO.write(records, args.output, 'fasta')



