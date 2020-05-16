#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Phylo
from calendar import isleap
import argparse
import generate_xml as gx
import subprocess
import shlex

parser = argparse.ArgumentParser()
# input files
parser.add_argument('--aln', help='input alignment file',
                    default='data/masked.fasta')
parser.add_argument('--initTree', help='input ML time tree',
                    default='data/tree.nwk')
parser.add_argument('--metadata', help='metadata file',
                    default='data/metadata.tsv')
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
                    default=(2020+((31)/366)))       # February 1
parser.add_argument('--s', help='number of susceptible individuals',
                    default=8883800)       # February 1
parser.add_argument('--latPeriod',
                    help='duration of exposed (latent) period (days)',
                    default=3)
parser.add_argument('--infPeriod', help='duration of infecitous period (days)',
                    default=5.5)
parser.add_argument('--P',
                    help='prop. of infection caused by I_high individuals',
                    default=0.8)
parser.add_argument('--R0M', help='desired mean for lognormal R0 prior',
                    default=1.5)
parser.add_argument('--R0S', help='desired s for lognormal R0 prior',
                    default=0.5)
parser.add_argument('--R0ChangeDate', help='date on which R0 changes by factor a',
                    default='NA')    # default is no R0 change
parser.add_argument('--EM',
                    default=1.0)
parser.add_argument('--IhM',
                    default=1.0E-8)
parser.add_argument('--exogInitM', help='mean on prior for inital size of exog',
                    default=1.0)
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
# Initial value ranges
parser.add_argument('--R0InitRange', help='range of initial R0 values',
                    default=(2, 5))
parser.add_argument('--exogGrInitRange', help='range for initial exog GR value',
                    default=(20, 30))
parser.add_argument('--exogInitRange', help='initial size of exog (at t0)',
                    default=(0, 1))
parser.add_argument('--importRateInitRange',
                    help='range of initial import rate', default=(5, 10))
parser.add_argument('--EInitRange', help='range of initial people in E class',
                    default=(1, 10))
# Run ID
parser.add_argument('--id', help='run id', default = 'Israel_SARS_CoV-2')
args = parser.parse_args()
args.base_path = args.aln.replace('.fasta', '')


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


def gen_xml(args):
    if args.R0ChangeDate != 'NA':
        name = str(round(args.sirT0, 3)) + '_' + str(args.P) + '_' +\
            str(round(args.R0ChangeDate, 3))
    else:
        name = str(round(args.sirT0, 3)) + '_' + str(args.P) + '_' +\
            args.R0ChangeDate
    args.out_name = args.id + name
    args.base_path = '/'.join(args.aln.split('/')[0:-1]) + '/fit_sirt0/'+args.id + \
        '_' + name
    gx.generate_xml(args)


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
    aln_name = args.aln.split('/')[-1].replace('.fasta', '')
    renamed_fasta_path = args.base_path+'/'+aln_name+'renamed.fasta'
    renamed_tree_path = args.base_path+'/'+aln_name+'renamed.newick'
    with open(renamed_fasta_path, 'w') as out_fasta:
        SeqIO.write(aligned, out_fasta, 'fasta')
    with open(renamed_tree_path, 'w') as out_tree:
        Phylo.write(tree, out_tree, 'newick')
    args.finalAln = renamed_fasta_path
    args.finalTree = renamed_tree_path
    return(args)


if __name__ == "__main__":
    args.importChangeDate = (2020+(31+29+20)/366)
    args.importChange = 0.25
    args.R0ChangeDate = (2020+(31+29+19)/366)
    args.xmlTemplate = 'config/SEIR_TEMPLATE.xml'
    for ph in [0.01, 0.02, 0.05, 0.07, 0.08, 0.09, 0.10]:
        args.pH = ph
        args.base_path = './data/beast/ph{0}'.format(str(ph))
        subprocess.run(shlex.split('mkdir -p {0}'.format(args.base_path)))
        args = process_fasta(args)
        args.out_name = 'ph{0}'.format(str(ph))
        gx.generate_xml(args)
