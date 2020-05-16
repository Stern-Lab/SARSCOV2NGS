#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Phylo
from calendar import isleap
import argparse
import generate_xml as gx

parser = argparse.ArgumentParser()
# input files
parser.add_argument('--aln', help='input alignment file',
                    default='data/mike/masked.fasta')
parser.add_argument('--initTree', help='input ML time tree',
                    default='data/mike/tree.nwk')
parser.add_argument('--metadata', help='metadata file',
                    default='data/mike/metadata.tsv')
parser.add_argument('--ref', help='reference fasta file',
                    default='data/MN908947.3.fasta')
parser.add_argument('--xmlTemplate', help='template xml file',
                    default='config/SEIR_TEMPLATE.xml')
# Input parameters (most of these are passed to PhyDyn XML generator)
parser.add_argument('--pH', help='proportion of I_high individuals',
                    default=0.2)
parser.add_argument('--T0', help='begining of time series',
                    default=2019.7)
parser.add_argument('--sirT0', help='begining of epi dynamics',
                    default=2020.085)       # February 1
parser.add_argument('--s', help='number of susceptible individuals',
                    default=9000000)       # February 1
parser.add_argument('--latPeriod',
                    help='duration of exposed (latent) period (days)',
                    default=4.2)
parser.add_argument('--infPeriod', help='duration of infecitous period (days)',
                    default=4.4)
parser.add_argument('--P',
                    help='prop. of infection caused by I_high individuals',
                    default=0.85)   # JLS superspreading paper, TODO find doi
parser.add_argument('--R0M', help='desired mean for lognormal R0 prior',
                    default=1.5)
parser.add_argument('--R0S', help='desired s for lognormal R0 prior',
                    default=0.5)
parser.add_argument('--R0Change', help='date on which R0 changes by factor a',
                    default='NA')    # default is no R0 change
parser.add_argument('--exogGrM', help='exog growth rate lognormal mean',
                    default=3.6)
parser.add_argument('--exogGrS',
                    help='exog growth rate lognormal standard deviation',
                    default=1)
parser.add_argument('--importRateM',
                    help='mean of import rate prior (exponential)', default=10)
parser.add_argument('--importRateU', help='upper limit of import rate',
                    default=10)
parser.add_argument('--importChangeDate',
                    help='date when importation changes', default=2050)
parser.add_argument('--importChange',
                    help='factor by which import rate changes', default=1)
parser.add_argument('--out_name', help='output file name',
                    default='Israel_SARS_CoV-2')
args = parser.parse_args()
args.base_dir = '/'.join(args.aln.split('/')[0:-1])
args.base_path = args.base_dir+'/'+args.out_name


# outputs date string for beast XML
def numeric_date(dt):
    days_in_year = 366 if isleap(dt.year) else 365
    res = dt.year + (dt.timetuple().tm_yday-0.5) / days_in_year
    return res


def rename_aln(aligned, pH, metadata):   # ph = probability of being in I_high
    names = {}
    for item in aligned:
        old_name = item.id
        names[old_name] = ''
        item_time = str(float(
            metadata.loc[metadata['strain'] == item.id]['numeric_date']))
        item.id = item.id.replace('_', '-')     # need to parse on '_' in beast
        item.id = item.id+'_'+str(item_time)
        if 'Israel' in item.id:
            if np.random.random(1)[0] > 1-pH:
                item.id = item.id + '_Ih'
            else:
                item.id = item.id + '_Il'
        else:
            item.id = item.id+'_exog'
        item.name = item.id
        item.description = item.id
        names[old_name] = item.id
    return(aligned, names)


def rescale_tree(tree):
    for non_terminal in tree.get_nonterminals():
        non_terminal.branch_length *= 1000
    for terminal in tree.get_terminals():
        terminal.branch_length *= 1000
    return(tree)


def process_fasta(args):
    metadata = pd.read_csv(args.metadata, sep='\t')
    # calcualte numeric date
    # removes rows with ambiguous dates
    exclude_dates = set(['2019', '2020', '2020-01', '2020-02', '2020-03',
                         '2020-01-XX', '2020-02-XX', '2020-03-XX',
                         '2020-04-XX'])
    metadata = metadata[~metadata['date'].isin(exclude_dates)]
    metadata['numeric_date'] = pd.to_datetime(
        metadata['date']).apply(numeric_date)
    aligned = list(SeqIO.parse(args.aln, "fasta"))
    aligned, names = rename_aln(aligned, args.pH, metadata)
    tree = Phylo.read(args.initTree, 'newick')
    tree = rescale_tree(tree)       # TODO WHY DO WE NEED THIS
    for tip in tree.get_terminals():
        tip.name = names[tip.name]
    renamed_fasta_path = args.base_path+'_renamed.fasta'
    renamed_tree_path = args.base_path+'_renamed.newick'
    with open(renamed_fasta_path, 'w') as out_fasta:
        SeqIO.write(aligned, out_fasta, 'fasta')
    with open(renamed_tree_path, 'w') as out_tree:
        Phylo.write(tree, out_tree, 'newick')
    args.finalAln = renamed_fasta_path
    args.finalTree = renamed_tree_path
    # TODO ADD TO DEFAULT PARAMETERS
    #TODO FIX DOIS IN TEMPLATE
    # First reported case on February 21
    # First peak in importation plot on Feb. 9
    args.sirT0 = (2020+((31+9)/366))
    # 362 reported cases outside china on Feb. 9
    args.exogInit = 10
    args.exogInitM = 20
    # Airports closed around 2020.189
    # Using March 20th (2020.219) based on importation analysis
    args.importChangeDate = (2020+(31+29+20)/366)
    # Should revisit this parameter
    args.importChange = 0.25
    # Constant R0
    args.base_path = args.base_path + '_constant_R0'
    args.out_name = args.out_name + '_constant_R0'
    gx.generate_xml(args)
    # R0 changes on date of school closures
    # Schools closure announced on March 12
    # Took effect on Sunday the 14th
    args.R0Change = (2020+(14+29+31)/366)
    args.base_path = args.base_path.replace('constant', 'piecewise')
    args.out_name = args.out_name.replace('constant', 'piecewise')
    gx.generate_xml(args)


if __name__ == "__main__":
    process_fasta(args)
