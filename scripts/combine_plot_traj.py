import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse
import subprocess
import shlex
import seaborn as sns
from calendar import isleap
import string
import baltic as bt


parser = argparse.ArgumentParser()
parser.add_argument('--log', default='Israel_SARS_CoV-2_constant_R0.log')
parser.add_argument('--traj', default='seir.Israel_SARS_CoV-2_constant_R0.traj')
parser.add_argument('--metadata', help='metadata file',
                    default='data/metadata.tsv')
parser.add_argument('--burnin', default=0.5)
parser.add_argument('--qs', help='quantiles for summary statistics', default='0.025,0.50,0.975')
args = parser.parse_args()
args.base_path = args.log.replace('.log', '')
args.qs = tuple([float(item) for item in args.qs.split(',')])

# Colors
polar_night = ["#2E3440", "#3b4252", "#434C5E", "#4C566A"]
dark_red = "#61272d"
country_col_dict = {'Africa': '#99A3B6',
                    'Asia': '#d070b9',
                    'Europe': '#d08770',
                    'Israel': '#bf616a',
                    'North America': '#8FBCBB',
                    'Oceania': '#5e81ac',
                    'South America': '#a3be8c'}

ph_col_dict = {0.01: '#81a1c1',
               0.02: '#ebcb8b',
               0.03: '#bf616a',
               0.04: '#a3be8c',
               0.04: '#5e81ac',
               0.05: '#bb8fbc',
               0.06: '#88c0d0',
               0.07: '#d08770',
               0.08: '#d070b9',
               0.09: '#BCBB8F',
               0.10: '#8FBCBB',
               0.20: '#be8ca3'}

# Matplotlib parameters
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
# Font colors
mpl.rcParams['text.color'] = polar_night[0]
mpl.rcParams['axes.labelcolor'] = polar_night[0]
mpl.rcParams['xtick.color'] = polar_night[0]
mpl.rcParams['ytick.color'] = polar_night[0]
# Font sizes
mpl.rcParams['figure.titlesize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['legend.title_fontsize'] = 16
# Border colors
mpl.rcParams['axes.edgecolor'] = polar_night[0]
# Changes seaborn violinplot edgecolor
sns.categorical._Old_Violin = sns.categorical._ViolinPlotter
class _My_ViolinPlotter(sns.categorical._Old_Violin):
    def __init__(self, *args, **kwargs):
        super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
        self.gray=polar_night[0]


sns.categorical._ViolinPlotter = _My_ViolinPlotter

# Combines beast logs
# not needed when only one replicate
def combine_runs(run_ids, args):
    # Trees
    combine_tree_str = ''
    combine_log_str = ''
    for run_id in run_ids:
        combine_tree_str += '-log {0}.trees '.format(run_id)
        combine_log_str += '-log {0}.log'.format(run_id)
    tree_logcombiner_cmd = 'logcombiner -b {0} -o {1} {2}'.format(
        args.burnin, args.base_path+'.tree', combine_tree_str)
    subprocess.run(shlex.split(tree_logcombiner_cmd))
    log_logcombiner_cmd = 'logcombiner -b {0} -o {1} {2}'.format(
        args.burnin, args.base_path+'.combined.log', combine_log_str)
    subprocess.run(shlex.split(log_logcombiner_cmd))
    args.burnin = 0
    return(args)


# Creates MCC tree from combined tree
def mcc_tree(tree_path, burnin):
    treeannotater_cmd = 'treeannotator -burnin {0} {1} {2}'.format(
        burnin, tree_path, tree_path.replace('.trees', '_mcc.tre'))
    subprocess.run(shlex.split(treeannotater_cmd))
    myTree = bt.loadNexus(tree_path.replace('.trees', '_mcc.tre'),
                          absoluteTime=False)
    # TODO hacky fix this later
    max_time = max([float(item.name.split('_')[1])
                    for item in myTree.Objects if item.branchType == 'leaf'])
    myTree.setAbsoluteTime(max_time)
    return(myTree)


# Plots tree
def plot_tree(tree, args):
    metadata = pd.read_csv(args.metadata, sep='\t')
    metadata['strain'] = metadata['strain'].str.replace('_', '-')
    def c_func(k):
        name = k.name.split('_')[0]
        item_metadata = metadata[metadata['strain'] == name]
        if item_metadata['country'].to_string(index=False).lstrip() == 'Israel':
            region = 'Israel'
        else:
            region = item_metadata['region'].to_string(index=False).lstrip()
        if region in country_col_dict.keys():
            return(country_col_dict[region])
        else:
            return(polar_night[-1])
    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8*2), constrained_layout=True)
    tree.plotTree(ax, x_attr=lambda k: k.absoluteTime,
                  colour_function=lambda k: polar_night[0],
                  branchWidth=lambda k: 1)
    tree.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
                    colour_function=lambda k: polar_night[0],
                    size_function=lambda k: 2*(50-30*k.height/tree.treeHeight))
    tree.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
                    colour_function=c_func,
                    size_function=lambda k: (50-30*k.height/tree.treeHeight))
    ax.set_xlim((2019.85, 2020.35))
    ax.set_yticks([])
    ax.ticklabel_format(useOffset=False)
    ax.grid(axis='x', ls='-', color='#d8dee9')
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    for i, region in enumerate(country_col_dict.keys()):
        ax.text(0.05, 0.5-i*0.035, region, color=country_col_dict[region],
                size=20, transform=ax.transAxes)
    fig.savefig('figures/mcc_tree_{0]}.pdf')
    plt.close(fig)


# Plots sample time histogram
def plot_sample_times(tree, args):
    times = pd.DataFrame([[item.name.split('_')[2],
                           item.name.split('_')[1]] for item in tree.Objects
                           if item.branchType=='leaf'])
    times.columns=['region', 'date']
    fig, axs = plt.subplots(1, 1, figsize=(6.4, 4.8), constrained_layout=True)
    sns.distplot(times[(times['region']=='Il') | (times['region']=='Ih')]['date'], color=country_col_dict['Israel'], ax=axs)
    sns.distplot(times[(times['region']=='exog')]['date'], color=country_col_dict['Africa'], ax=axs)
    axs.text(0.05, 0.9, 'Israel', color=country_col_dict['Israel'],
        size=20, transform=axs.transAxes)
    axs.text(0.05, 0.9-0.07, 'Global', color=country_col_dict['Africa'],
        size=20, transform=axs.transAxes)
    axs.set_xlabel('Sampling date (year)')
    axs.set_xlim((2019.85, 2020.35))
    axs.ticklabel_format(useOffset=False)
    axs.set_ylabel('Count')
    fig.savefig(args.log.replace('.log', '_dates.pdf'))
    plt.close(fig)


# Turns pandas datetime object into numeric date
def numeric_date(dt):
    days_in_year = 366 if isleap(dt.year) else 365
    res = dt.year + (dt.timetuple().tm_yday-0.5) / days_in_year
    return res


# Calculates quantiles of time series
def process_traj(traj, item, args):
    item_processed = traj[['t', item]].groupby('t').quantile(args.qs).unstack()
    item_processed.columns = args.qs
    return(item_processed)


# formats results for output CSV
def add_formatted_result(formatted_results, key, item):
    new_result = [key,
                  '%s' % float('%.5g' % item.iloc[1]) + ' (' +
                  '%s' % float('%.5g' % item.iloc[0]) + ', ' +
                  '%s' % float('%.5g' % item.iloc[2]) + ')']
    formatted_results.append(new_result)
    return(formatted_results)


# reads in beast logs and phydyn trajectory files, calculates quantiles,
# creates formatted output csv
def process_results(args):
    results = {}
    formatted_results = []
    log = pd.read_csv(args.log, sep="\t", comment="#")
    log = log[math.ceil(len(log)*args.burnin):]
    for d in log.columns[1:]:   # excludes sample column
        results[d] = log[d].quantile(args.qs)
    results['seir.R0.all'] = log['seir.R0']
    results['seir.a.all'] = log['seir.a']
    for key, item in results.items():
        formatted_results = add_formatted_result(formatted_results, key, item)
    # reads in PhyDyn trajectories
    traj = pd.read_csv(args.traj, sep='\t')
    for item in traj.columns[2:]:
        results[item] = process_traj(traj, item, args)
    max_time = max(traj['t'])
    # endpoint # of recovered individuals
    results['tot_R.all'] = traj[traj['t'] == max(traj['t'])]['R']
    # endpoint cumulative # of infectious individuals
    if 'infectious' in traj.columns:
        results['tot_infectious.all'] = traj[
            traj['t'] == max(traj['t'])]['infectious']
    # I includes all INFECTIOUS individuals (E+Il+Ih)
    results['I'] = results['Il'].add(results['Ih'], fill_value=0)
    formatted_results = add_formatted_result(formatted_results,
                                             'cumulative incidence',
                                             results['infections'].iloc[-1, :])
    formatted_results = pd.DataFrame(formatted_results,
                                     columns=['Parameter', 'Median (95% HPD)'])
    formatted_results.to_csv(args.log.replace('.log', '.csv'), index=False)
    return(results, formatted_results)


# Individual trajectory plots for each p_h value
def plot_traj_individual(all_results, israel_reported_data, plot_ph):
    fig, axs = plt.subplots(len(all_results), 2,
                            figsize=(6.4*2, 4.8*len(all_results)),
                            constrained_layout=True)
    for i, ph in enumerate(plot_ph):
        axs[i][0].plot(all_results[ph]['I'].index,
                       all_results[ph]['I'][0.5],
                       linewidth=3,
                       color=ph_col_dict[ph])
        axs[i][0].fill_between(all_results[ph]['I'].index,
                               all_results[ph]['I'][0.025],
                               all_results[ph]['I'][0.975],
                               color=ph_col_dict[ph],
                               alpha=0.25)
        axs[i][0].set_ylabel('Infectious individuals')
        axs[i][0].text(0.05, 0.9-0.07, '$p_h$ = {0}'.format(ph),
                       color=ph_col_dict[ph],
                       size=20, transform=axs[i][0].transAxes)
        sns.scatterplot(x="numeric_date", y="cumulative cases",
                        data=israel_reported_data,
                        color=polar_night[0], ax=axs[i][1])
        axs[i][1].plot(all_results[ph]['R'].index,
                       all_results[ph]['R'][0.5],
                       linewidth=3,
                       color=ph_col_dict[ph])
        axs[i][1].fill_between(all_results[ph]['R'].index,
                               all_results[ph]['R'][0.025],
                               all_results[ph]['R'][0.975],
                               color=ph_col_dict[ph],
                               alpha=0.25)
        axs[i][1].set_ylabel('Cumulative incidence')
        for j, ax in enumerate(axs[i]):
            ax.ticklabel_format(useOffset=False)
            ax.set_yscale('log')
            ax.set_xlabel('Month (2020)')
            ax.set_ylim(1, )
            ax.set_xlim(2020+(31+1)/366-0.035,
                        2020+(31+29+31+31+1)/366+0.035)
            ax.set_xticks((2020+(31+1)/366, 2020+(31+29+1)/366,
                           2020+(31+29+31+1)/366, 2020+(31+29+31+31+1)/366))
            ax.set_xticklabels(('Feb.', 'March', 'April', 'May'))
            ax.text(-0.15, 1, string.ascii_uppercase[i*2 + j],
                    color=polar_night[0], size=24,
                    transform=ax.transAxes, fontweight='bold')
        fig.savefig('figures/individual_trajectories.pdf')
        plt.close(fig)


# Plots trajectory of selected p_h values on the same plot
def plot_traj_combined(all_results, israel_reported_data, plot_ph):
    fig, axs = plt.subplots(1, 2, figsize=(6.4*2, 4.8),
                            constrained_layout=True)
    sns.scatterplot(x="numeric_date", y="cumulative cases",
                    data=israel_reported_data,
                    color=polar_night[0], ax=axs[1])
    for i, ph in enumerate(plot_ph):
        axs[0].plot(all_results[ph]['I'].index,
                    all_results[ph]['I'][0.5],
                    linewidth=3,
                    color=ph_col_dict[ph])
        axs[0].fill_between(all_results[ph]['I'].index,
                            all_results[ph]['I'][0.025],
                            all_results[ph]['I'][0.975],
                            color=ph_col_dict[ph],
                            alpha=0.25)
        axs[0].set_ylabel('Infectious individuals')
        axs[0].text(0.05, 0.9-i*0.07, '$p_h$ = {0}'.format(ph),
                    color=ph_col_dict[ph],
                    size=20, transform=axs[0].transAxes)
        axs[1].plot(all_results[ph]['R'].index,
                    all_results[ph]['R'][0.5],
                    linewidth=3,
                    color=ph_col_dict[ph])
        axs[1].fill_between(all_results[ph]['R'].index,
                            all_results[ph]['R'][0.025],
                            all_results[ph]['R'][0.975],
                            color=ph_col_dict[ph],
                            alpha=0.25)
        axs[1].set_ylabel('Cumulative incidence')
    for i, ax in enumerate(axs):
        ax.ticklabel_format(useOffset=False)
        ax.set_yscale('log')
        ax.set_xlabel('Month (2020)')
        ax.set_ylim(1, )
        ax.set_xlim(2020+(31+1)/366-0.035,
                    2020+(31+29+31+31+1)/366+0.035)
        ax.set_xticks((2020+(31+1)/366, 2020+(31+29+1)/366,
                       2020+(31+29+31+1)/366, 2020+(31+29+31+31+1)/366))
        ax.set_xticklabels(('Feb.', 'March', 'April', 'May'))
        ax.text(-0.15, 1, string.ascii_uppercase[i], color=polar_night[0],
                size=24, transform=ax.transAxes, fontweight='bold')
    fig.savefig('figures/combined_trajectory.pdf')
    plt.close(fig)


# Plots epi parameters across a range of p_h values
def plot_ph_sensitivity(R0_s, a_s, R_s, israel_reported_data):
    colors = [item for item in ph_col_dict.values()]
    sns.set_palette(sns.color_palette(colors))
    fig, axs = plt.subplots(1, 3, figsize=(6.4*3, 4.8),
                            constrained_layout=True)
    sns.violinplot(x='ph', y='seir.R0', data=R0_s, ax=axs[0], saturation=1)
    axs[0].set_ylabel("$R_0$")
    axs[0].set_ylim(0, 3)
    sns.violinplot(x='ph', y='seir.a', data=a_s, ax=axs[1], saturation=1)
    axs[1].set_ylabel(r"$ \alpha $")
    axs[1].set_xlabel("$p_h$")
    axs[1].set_ylim(0, 1)
    sns.violinplot(x='ph', y='R', data=R_s, ax=axs[2], saturation=1)
    axs[2].axhline(y=max(israel_reported_data['cumulative cases']),
                   linestyle='--', color=polar_night[0])
    axs[2].set_ylabel('Cumulative infectious')
    axs[2].set_yscale('log')
    axs[2].set_ylim(1, )
    for i, ax in enumerate(axs):
        ax.set_xlabel("$p_h$")
        ax.text(-0.15, 1, string.ascii_uppercase[i], color=polar_night[0],
                size=24, transform=ax.transAxes, fontweight='bold')
    fig.savefig('figures/ph_sensitivity.pdf')
    plt.close(fig)


def process_reported_data(data_path, all_results):
    reported_data = pd.read_csv(data_path)
    israel_reported_data = reported_data[
        reported_data['countriesAndTerritories'] == 'Israel']
    israel_reported_data['numeric_date'] = pd.to_datetime(
        israel_reported_data['dateRep'], dayfirst=True).apply(numeric_date)
    israel_reported_data.sort_values(by='numeric_date', inplace=True)
    israel_reported_data['cumulative cases'] = israel_reported_data[
        'cases'].cumsum()
    israel_reported_data = israel_reported_data[israel_reported_data[
        'numeric_date'] <= max(all_results[0.01]['I'].index)]
    return(israel_reported_data)


if __name__ == "__main__":
    subprocess.run(shlex.split('mkdir figures'))
    all_results = {}
    R0_s = pd.DataFrame()
    a_s = pd.DataFrame()
    R_s = pd.DataFrame()
    ph_s = [0.01, 0.02, 0.05, 0.07, 0.08, 0.09, 0.10]
    for ph in ph_s:
        args.log = 'data/beast/ph{0}/ph{0}.log'.format(str(ph))
        args.traj = 'data/beast/ph{0}/seir.ph{0}.traj'.format(str(ph))
        tree_path = 'data/beast/ph{0}/ph{0}.trees'.format(str(ph))
        myTree = mcc_tree(tree_path, args.burnin*100)
        plot_tree(myTree, args)
        plot_sample_times(myTree, args)
        results, formatted_results = process_results(args)
        all_results[ph] = results
        R0 = pd.DataFrame(results['seir.R0.all'])
        # removes values outside of 95% HPD
        R0 = R0[R0['seir.R0'].between(R0.quantile(args.qs[0])[0],
                                      R0.quantile(args.qs[1])[0])]
        R0['ph'] = ph
        R0_s = R0_s.append(R0)
        a = pd.DataFrame(results['seir.a.all'])
        a['ph'] = ph
        a = a[a['seir.a'].between(a.quantile(args.qs[0])[0],
                                  a.quantile(args.qs[2])[0])]
        a_s = a_s.append(a)
        R = pd.DataFrame(results['tot_R.all'])
        R = R[R['R'].between(R.quantile(args.qs[0])[0],
                             R.quantile(args.qs[2])[0])]
        R['ph'] = ph
        R_s = R_s.append(R)
    israel_reported_data = process_reported_data('data/ecdc_reported.csv',
                                                 all_results)
    plot_ph_sensitivity(R0_s, a_s, R_s, israel_reported_data)
    plot_traj_individual(all_results, israel_reported_data, ph_s)
    plot_ph = [0.01, 0.05, 0.10]
    plot_traj_combined(all_results, israel_reported_data, plot_ph)
