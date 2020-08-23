#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Phylo
from calendar import isleap
import argparse
from scripts import generate_xml as gx
import subprocess
import shlex
from scipy import integrate


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


# rescales branch lengths on tree to years
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
    aligned, names = rename_aln(aligned, args.ph, metadata)
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


def main():
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
    parser.add_argument('--ph', help='proportion of I_high individuals',
                        default=0.2)
    parser.add_argument('--T0', help='begining of time series',
                        default=2019.7)
    parser.add_argument('--sirT0', help='begining of epi dynamics',
                        default=2020+(31+1)/366)       # February 1
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
                        default=(2020+(31+29+19)/366))  
    parser.add_argument('--EM', help='mean of initial number of people in E',
                        default=1.0)
    parser.add_argument('--exogInitM', help='mean on prior for inital size of exog',
                        default=1.0)
    parser.add_argument('--exogGR', help='growth rate of exog', default=24.0)
    parser.add_argument('--importRate',
                        help='import rate', default=1000)
    parser.add_argument('--id', help='run id', default = 'Israel_SARS_CoV-2')
    args = parser.parse_args()
    # This generates the XMLs used for initial journal submission
    run_id = 'june3_fix_import'
    args.xmlTemplate = 'config/SEIR_TEMPLATE_FIX_IMPORT.xml'
    subprocess.run(shlex.split('mkdir {0}'.format(run_id)))
    for ph in [0.02, 0.05, 0.1, 0.2, 0.5, 0.80]:
        for eta in [10, 100, 1000, 2500, 5000]:
            args.importRate = eta
            args.xmlTemplate = 'config/SEIR_TEMPLATE_FIX_IMPORT.xml'
            args.base_path = '{0}/{0}_{1}_{2}'.format(run_id, ph, eta)
            args.out_name = '{0}_{1}_{2}'.format(run_id, ph, eta)
            subprocess.run(shlex.split('mkdir {0}'.format(args.base_path)))
            args = process_fasta(args)
            gx.generate_xml(args)
    # This generates the XMLs used in the updated analysis with fixed eta
    # This XML template uses coupled MCMC to help convergence
    # Also adds an up-down operator on R0/E
    run_id = 'june11_CMCMC_Final'
    args.xmlTemplate = 'config/SEIR_TEMPLATE_FIX_IMPORT_UPDOWN_CMCMC'
    for ph in [0.02, 0.05, 0.1, 0.2, 0.5, 0.80]:
        for eta in [10, 100, 1000, 2500, 5000]:
            args.importRate = eta
            args.base_path = '{0}/{0}_{1}_{2}'.format(run_id, ph, eta)
            args.out_name = '{0}_{1}_{2}'.format(run_id, ph, eta)
            subprocess.run(shlex.split('mkdir {0}'.format(args.base_path)))
            args = process_fasta(args)
            gx.generate_xml(args)
    # Finally, this generates the XMLs in the updated analysis using an empirical curve
    # for eta. This is estimated based on a global tree. A piecewise exponential
    # function is fit to the TMRCA estimates into Israel using the scripts/import_analysis.py file
    # still uses the up-down operator on R0/E and uses coupled MCMC
    run_id = 'july22_TMRCA_Import_Final'
    subprocess.run(shlex.split('mkdir data/{0}'.format(run_id)))
    args.xmlTemplate = 'config/SEIR_TEMPLATE_IMPORTDEF.xml'
    import_rates = pd.read_csv('data/import_rates.csv')
    import_rates = {item['theta']: item['eta'] for i, item in import_rates.iterrows()}
        for ph in [0.02, 0.1]:
            for theta in [0.8, 0.9]:
                args.ph = ph
                args.importRate = import_rates[theta]
                args.base_path = 'data/{0}/{0}_{1}_{2}'.format(run_id, ph, theta)
                args.out_name = '{0}_{1}_{2}'.format(run_id, ph, theta)
                subprocess.run(shlex.split('mkdir {0}'.format(args.base_path)))
                args = process_fasta(args)
                gx.generate_xml(args)


if __name__ == "__main__":
    main()