import argparse

import PyAstronomy
import numpy as np
from Bio import Phylo
from tqdm import tqdm
import os
from ete3 import Tree
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from process_beast_log import numericToCalander

def numericToCalander(d):
    return PyAstronomy.pyasl.decimalYearGregorianDate(d)

def get_importation_dates(treepath, branch_dates):
    """
    get importation dates according to Volz report here:
    http://sarscov2phylodynamics.org/2020/04/08/importations.html
    :param treepath: a path to a tree file
    :param branch_dates: json file containing dates for each internal node.
    one of NextStrain standard outputs
    :return: a list containing the distribution of dates.
    """
    # open tree files
    tree = Phylo.read(treepath, "newick")
    with open(branch_dates) as json_file:
        nodes_to_dates = json.load(json_file)

    tips = [t for t in tree.get_terminals() if 'Israel' in t.name]
    dates_distribution = []
    for t in tips:
        tip_date = nodes_to_dates['nodes'][t.name]['numdate']
        node_path = [n for n in tree.get_path(t) if 'Wuhan' not in n.name or 'Israel' not in n.name]
        # remove outlier nodes
        #node_path = [n for n in node_path if (n.name != 'NODE_0000034' or n.name != 'NODE_0000035')]
        node_path = [n for n in node_path]
        node_names = [n.name for n in node_path if '.' not in n.name] # remove the 1.00 addition after polytomy splitting
        #dates = [nodes_to_dates['nodes'][n]['numdate'] for n in node_names]
        dates = []
        for node in node_names:
            if nodes_to_dates['nodes'][node]['numdate'] > 2019.9:
                print("WARNINIG: Detected a date beroe 2019.9")
            dates.append(nodes_to_dates['nodes'][node]['numdate'])
        mean_dates = [np.mean([tip_date, node_date]) for node_date in dates]
        dates_distribution.extend(mean_dates)
    return dates_distribution

def analyze_simulated_trees():

    n = 20
    dfs = []
    for i in range(n):
        treetime_outfile = fr"/Users/daniellemiller/Documents/GitHub/COVID19/data/" \
                           fr"importations_sampled_trees/sim_treetime_{i+1}.nwk"
        node_data_outfile = fr"/Users/daniellemiller/Documents/GitHub/COVID19/data/" \
                            fr"importations_sampled_trees/sim_branch_lengths_{i+1}.json"
        vals = get_importation_dates(treetime_outfile, node_data_outfile)
        df = pd.DataFrame({'sim': i+1, 'dates': vals})
        dfs.append(df)
    res = pd.concat(dfs)

    # plot the kde
    ax = sns.distplot(res['dates'], hist=False, kde_kws={'shade': True}, color='#D9C12B')
    ax.axvline(2020.14, label='2020-02-21', linestyle='--', color='red', alpha=0.5)
    ax.axvline(2020.253, label='last sample', linestyle='--')
    # rename the xticks with labels
    x_ticks = ax.get_xticks()
    ax.set_xticks(x_ticks[::2])
    xlabels = [pd.to_datetime(numericToCalander(y)).date() for y in x_ticks[::2]]
    ax.set_xticklabels(xlabels)
    plt.legend()
    sns.despine(offset=15)
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree_input", required=True, help="The raw ML tree")
    parser.add_argument("--aln", required=True, help="alignment file path")
    parser.add_argument("--metadata", required=True, help="metadata tsv file path")
    parser.add_argument("--nsim",type=int, default=20, help="number of simulations")
    args = parser.parse_args()

    skeleton_tree_path = args.tree_input
    aln = args.aln
    metadata = args.metadata
    nsim = args.nsim

    for i in tqdm(range(nsim)):
        out_file = fr"../results/sim_tree_{i+1}.nwk"
        # resolve polytomies and collapse low confidence splits

        tree = Tree(skeleton_tree_path)
        tree.resolve_polytomy(recursive=True)
        tree.write(outfile=out_file)


        # re-run tree time to create a new pull of trees
        treetime_outfile = fr"../results/sim_treetime_{i+1}.nwk"
        node_data_outfile = fr"../results/sim_branch_lengths_{i+1}.json"
        clock_rate = np.random.uniform(0.0009,0.0015)
        os.system(f"augur refine \
                        --tree {out_file} \
                        --alignment {aln} \
                        --metadata {metadata} \
                        --output-tree {treetime_outfile} \
                        --output-node-data {node_data_outfile} \
                        --root Wuhan-Hu-1/2019 \
                        --timetree \
                        --clock-rate {clock_rate} \
                        --coalescent skyline \
                        --date-inference marginal \
                        --divergence-unit mutations \
                        --date-confidence \
                        --no-covariance")
        #re-run ancestral state reconstruction
        ancestral_node_data_outfile = fr"/Users/daniellemiller/Documents/GitHub/COVID19/data/" \
                           fr"importations_sampled_trees/sim_nt_muts_{i+1}.json"
        os.system(f"augur ancestral \
                        --tree {treetime_outfile} \
                        --alignment {aln} \
                        --output-node-data {ancestral_node_data_outfile} \
                        --inference joint \
                        --keep-ambiguous")





