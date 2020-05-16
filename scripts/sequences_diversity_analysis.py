import argparse

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, Phylo
import itertools


def hamming_pairwise_distance(sequences, country='Israel'):
    """
    path to a fasta file containing sequences, pairwise distance will be calculated by contry
    :return: a dictionary mapping records to hamming distances
    """
    records = []
    for rec in SeqIO.parse(sequences,'fasta'):
        if country.lower() in rec.id.lower():
            records.append(rec)
    # get all pairs combinations by ids
    hamming_dict = {}
    pairs = list(itertools.combinations([r.id for r in records], 2))
    for rec1, rec2 in pairs:
        s1 = str([t for t in records if t.id == rec1][0].seq).lower()
        s2 = str([t for t in records if t.id == rec2][0].seq).lower()
        distance = sum(c1 != c2 for c1, c2 in zip(s1, s2) if c1 != 'n' and c2 != 'n')
        key = tuple(sorted([rec1,rec2]))    # sort key dictionary for an easy mapping of ids to distance
        hamming_dict[key] = distance
    return hamming_dict

def pairwise_delta_t(sequences, metadata, country='Israel'):
    """
    for each two sequences from a country return their difference in sampling dates in days
    :param sequences: a file containing sequences
    :param metadata: a tsv file holding dates for each sequence in sequences
    :param country: string representing the country of interest
    :return: a dictionary mapping records to date differences
    """
    records = []
    for rec in SeqIO.parse(sequences,'fasta'):
        if country.lower() in rec.id.lower():
            records.append(rec.id)
    df = pd.read_table(metadata)
    df = df[df['strain'].isin(records)] # filter the data frame to contain only relevant sequences - huge df originally
    delta_t_dict = {}

    pairs = list(itertools.combinations(records, 2))
    for p1, p2 in pairs:
        date1 = df[df['strain'] == p1]['date'].values[0]
        date2 = df[df['strain'] == p2]['date'].values[0]
        delta = abs((pd.to_datetime(date1, format="%d/%m/%Y") - pd.to_datetime(date2, format="%d/%m/%Y")).days)
        key = tuple(sorted([p1, p2]))
        delta_t_dict[key] = delta
    return delta_t_dict

def branch_length_pairwise_distance(sequences, tree, country='Israel'):
    """
    path to a fasta file containing sequences, pairwise branch length distance will be calculated by contry
    :return: a dictionary mapping records to hamming distances
    """
    records = []
    for rec in SeqIO.parse(sequences,'fasta'):
        if country.lower() in rec.id.lower():
            records.append(rec.id)
    # get all pairs combinations by ids
    bl_dict = {}
    pairs = list(itertools.combinations(records, 2))

    # read raw tree file
    tree = Phylo.read(tree, 'newick')
    tips = [t for t in tree.get_terminals() if country.lower() in t.name.lower()]

    for rec1, rec2 in pairs:
        tip1 = [t for t in tips if t.name == rec1][0]
        tip2 = [t for t in tips if t.name == rec2][0]
        distance = tree.distance(tip1, tip2)
        key = tuple(sorted([rec1,rec2]))    # sort key dictionary for an easy mapping of ids to distance
        bl_dict[key] = distance
    return bl_dict

def plot_regression(recs_2_days, recs_2_distance, distance_metric='hamming', out=None):
    """
    plot a regression plot
    :param recs_2_days:
    :param recs_2_distance:
    :param distance_metric:
    :return:
    """
    time_and_distance=[]
    for key in recs_2_days:
        time_and_distance.append((recs_2_days[key], recs_2_distance[key]))
    df = pd.DataFrame(time_and_distance,columns=['time', 'distance'])

    g = df.groupby('time').median().reset_index()
    g2 = df.groupby('time').size().reset_index(name='size')
    merged = pd.merge(g, g2, on='time')

    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        sns.regplot(x='time', y='distance', data=merged, color='#E38074', line_kws={'color': 'skyblue'},
                scatter_kws={'s': merged['size']})
        plt.xlabel(r'$\Delta$ days between samples')
        if distance_metric == 'hamming':
            plt.ylabel("Pairwise Hamming distance")
        else:
            plt.ylabel("Pairwise Branch length distance")
        plt.tight_layout()
        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.show()


def plot_boxplot(recs_2_days, recs_2_distance, distance_metric='hamming', out=None):
    """
    plot a regression plot
    :param recs_2_days:
    :param recs_2_distance:
    :param distance_metric:
    :return:
    """
    time_and_distance=[]
    for key in recs_2_days:
        time_and_distance.append((recs_2_days[key], recs_2_distance[key]))
    df = pd.DataFrame(time_and_distance,columns=['time', 'distance'])

    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        sns.boxplot(x='time', y='distance', data=df, color = 'white', fliersize = 0)
        plt.xlabel(r'$\Delta$ days between samples')
        if distance_metric == 'hamming':
            plt.ylabel("Pairwise Hamming distance")
        else:
            plt.ylabel("Pairwise Branch length distance")
        plt.tight_layout()
        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="The raw ML tree")
    parser.add_argument("--aln", required=True, help="alignment file path")
    parser.add_argument("--metadata", required=True, help="metadata tsv file path")
    parser.add_argument("--country", default="Israel", help="country")
    args = parser.parse_args()

    sequences = args.aln
    metadata = args.metadata
    tree = args.tree

    country = args.country
    delta_days = pairwise_delta_t(sequences, metadata, country)
    pairwise_hamming = hamming_pairwise_distance(sequences, country)
    pairwise_bl = branch_length_pairwise_distance(sequences, tree, country)

    plot_regression(delta_days, pairwise_hamming, distance_metric='hamming', out='../results/days_2_hamming.png')
    plot_regression(delta_days, pairwise_bl, distance_metric='BL', out='../results/days_2_BL.png')

    plot_boxplot(delta_days, pairwise_hamming, distance_metric='hamming', out='../results/box_days_2_hamming.png')
    plot_boxplot(delta_days, pairwise_bl, distance_metric='BL', out='../results/box_days_2_BL.png')











