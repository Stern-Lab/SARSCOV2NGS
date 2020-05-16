import argparse
from Bio import SeqIO, Phylo
import pandas as pd
import numpy as np
from tqdm import tqdm


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



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences after filtration")
    parser.add_argument("--exogenous", type = int, default=105, help="number of exogenous sequences to sample")
    parser.add_argument("--metadata", type=str, required=True, help="path to a metadata file")
    parser.add_argument("--tree", type=str, required=True, help="path to a tree file")
    parser.add_argument("--country", type=str, default='Israel', help="country identifier")
    parser.add_argument("--output", required=True, help="FASTA output file")
    args = parser.parse_args()

    metadata = pd.read_table(args.metadata)
    tree = Phylo.read(args.tree, "newick")
    # make sure the tree is rooted
    tree.root_with_outgroup({'name': 'Wuhan-Hu-1/2019'})
    country = args.country.lower()
    metadata['country'] = metadata['country'].apply(lambda x: x.lower())

    sampled_strains = time_sampling(args.sequences, metadata, country, args.exogenous)

    # # add to the list the strains from closest match
    # with open(r'/Volumes/STERNADILABHOME$/volume3/COVID19/data/sampling/exclude.txt', 'r') as f:
    #     exclude = [s.strip() for s in f.readlines()]
    # sampled_strains = [s for s in sampled_strains if s not in exclude]

    sampled = list(set(sampled_strains + ['Wuhan-Hu-1/2019']))    # add the root

    records = [rec for rec in SeqIO.parse(args.sequences, 'fasta') if rec.id in sampled]
    cnt = SeqIO.write(records, args.output, 'fasta')





