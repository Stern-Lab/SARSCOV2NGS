import argparse
import numpy as np
import os
from Bio import SeqIO
import re
import scipy.stats


def world_bootstrap(sequences, outpath, country="israel"):
    """
    sample with replacement from an alignment
    """
    records = [rec for rec in SeqIO.parse(sequences, "fasta")]
    wuhan = [rec for rec in records if rec.name == 'Wuhan-Hu-1/2019']
    world_records = [rec for rec in records if country not in rec.name.lower()]
    country_records = [rec for rec in records if country in rec.name.lower()]
    n = len(world_records)
    sampled_idx = list(np.random.choice(n, n))  # running times
    sampled = [world_records[i] for i in set(sampled_idx)]

    recs_2_write = sampled + country_records + wuhan
    cnt = SeqIO.write(recs_2_write, outpath, "fasta")
    return cnt


def country_bootstrap(sequences, outpath, country="israel"):
    """
    sample with replacement from an alignment only israeli sequences
    """
    records = [rec for rec in SeqIO.parse(sequences, "fasta")]
    wuhan = [rec for rec in records if rec.name == 'Wuhan-Hu-1/2019']
    world_records = [rec for rec in records if country not in rec.name.lower()]
    country_records = [rec for rec in records if country in rec.name.lower()]
    n = len(country_records)
    sampled_idx = list(np.random.choice(n, n))  # running times
    sampled = [country_records[i] for i in set(sampled_idx)]

    recs_2_write = sampled + world_records + wuhan
    cnt = SeqIO.write(recs_2_write, outpath, "fasta")
    return cnt


def write_snakefile(filepath, base_snakefile):
    """
    create a snake file on each directory based on a template
    """
    with open(base_snakefile, 'r') as base_file:
        text = base_file.read()
        text_2_write = re.sub(r'(FOLDER\n)',f'FOLDER = "world/{os.path.basename(filepath)}/"\n', text)
    with open(os.path.join(filepath,"Snakefile"), 'w') as out_file:
        cnt = out_file.write(text_2_write)
    return cnt


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Down sample sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--type", required=True, help="world \ country")
    parser.add_argument("--aln", required=True, help="alignment file path")
    parser.add_argument("--working_dir", required=True, help="working dir file path")
    parser.add_argument("--n", type=int, default=1000, help="number of bootstrap replicates")
    args = parser.parse_args()

    base_dir = args.working_dir
    aln = args.aln
    num_samples = args.n
    i = args.i

    # updated version on cluster - distributed from each index (can ajust to a for loop - more time consuming)
    alias = args.type
    if os.path.exists(os.path.join(base_dir, f"{alias}_{i}_bootstrap")):
        if os.path.exists(os.path.join(base_dir, f"{alias}_{i}_bootstrap/data")):
            outpath = os.path.join(base_dir, f"{alias}_{i}_bootstrap/data/sequences.fasta")
            if "world" in alias:
                world_bootstrap(aln, outpath)
            else:
                country_bootstrap()
    else:
        os.makedirs(os.path.join(base_dir, f"{alias}_{i}_bootstrap"))
        os.makedirs(os.path.join(base_dir, f"{alias}_{i}_bootstrap/data"))
        outpath = os.path.join(base_dir, f"{alias}_{i}_bootstrap/data/sequences.fasta")
        if "world" in alias:
            world_bootstrap(aln, outpath)
        else:
            country_bootstrap()





