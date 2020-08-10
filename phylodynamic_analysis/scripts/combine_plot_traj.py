#!/usr/bin/env python3
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import numpy as np
import argparse
import subprocess
import shlex
import seaborn as sns
from calendar import isleap
import string
from scripts import sars_cov2_plotting_style as plot_style
import baltic as bt
import pickle
import csv


# Creates MCC tree from combined tree
def calc_mcc_tree(tree_path, burnin):
    treeannotater_cmd = 'treeannotator -burnin {0} {1} {2}'.format(
        burnin, tree_path, tree_path.replace('.trees', '_mcc.tre'))
    subprocess.run(shlex.split(treeannotater_cmd))
    myTree = bt.loadNexus(tree_path.replace('.trees', '_mcc.tre'),
                          absoluteTime=False)
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
        if region in plot_style.country_col_dict.keys():
            return(plot_style.country_col_dict[region])
        else:
            return(plot_style.polar_night[-1])
    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8*2), constrained_layout=True)
    tree.plotTree(ax, x_attr=lambda k: k.absoluteTime,
                  colour_function=lambda k: plot_style.polar_night[0],
                  branchWidth=lambda k: 1)
    tree.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
                    colour_function=lambda k: plot_style.polar_night[0],
                    size_function=lambda k: 2*(50-30*k.height/tree.treeHeight))
    tree.plotPoints(ax, x_attr=lambda k: k.absoluteTime,
                    colour_function=c_func,
                    size_function=lambda k: (50-30*k.height/tree.treeHeight))
    ax.set_xlim((2019.85, 2020.35))
    ax.set_xticks(pd.to_datetime(['12-01-19', '01-01-20', '02-01-20', '03-01-20', '04-01-20']).\
                    to_series().apply(numeric_date).to_list())
    ax.set_xticklabels(['Dec.\n2019', 'Jan.\n2020', 'Feb.\n2020', 'March\n2020', 'April\n2020'])
    ax.set_yticks([])
    ax.grid(axis='x', ls='-', color='#d8dee9')
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    for i, region in enumerate(plot_style.country_col_dict.keys()):
        ax.text(0.05, 0.5-i*0.035, region, color=plot_style.country_col_dict[region],
                size=20, transform=ax.transAxes)
    fig.savefig(f'figures/{args.plot_name}_{args.ph}_{args.eta}_mcc_tre.pdf')
    plt.close(fig)


# Plots sample time histogram
def plot_sample_times(tree, args):
    times = pd.DataFrame([[item.name.split('_')[2],
                           item.name.split('_')[1]] for item in tree.Objects
                           if item.branchType=='leaf'])
    times.columns=['region', 'date']
    fig, axs = plt.subplots(1, 1, figsize=(6.4, 4.8), constrained_layout=True)
    sns.distplot(times[(times['region']=='Il') | (times['region']=='Ih')]['date'], color=plot_style.country_col_dict['Israel'], ax=axs)
    sns.distplot(times[(times['region']=='exog')]['date'], color=plot_style.country_col_dict['Africa'], ax=axs)
    axs.text(0.05, 0.9, 'Israel', color=plot_style.country_col_dict['Israel'],
        size=20, transform=axs.transAxes)
    axs.text(0.05, 0.9-0.07, 'Global', color=plot_style.country_col_dict['Africa'],
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
def process_traj(results, max_time, traj, item, args):
    hpd_width = args.qs[2] - args.qs[0]
    item_traj = traj[['t', item]]
    item_cum = item_traj[item_traj['t'] == max(item_traj['t'])]
    cum_hpd_limits = hpd(item_cum[item], args.qs)
    cum_hpd = item_cum[item_cum[item].between(cum_hpd_limits[0],
                                              cum_hpd_limits[2])]
    results[item+'_cum_HPD'] = cum_hpd
    item_traj_hpd = pd.DataFrame(item_traj.groupby('t')[item].apply(hpd, args.qs).to_list(),
                                 columns=args.qs)
    item_traj_hpd.index = item_traj['t'][0:len(item_traj_hpd)]
    results[item+'_traj_HPD'] = item_traj_hpd
    return(results)


# mostly from arviz
# but works directly with a pandas object
# and also returns median
# https://github.com/arviz-devs/arviz/blob/master/arviz/stats/stats.py
def hpd(dat, qs):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = dat.sort_values()
    # length of data
    n = len(dat)
    # number of values we are keeping
    # thus, the index of the begining of 
    # the HPD interval must be <= this value
    # this gives us the tail of the distribution
    interval_idx_inc = int(np.floor(width * n))
    # number of values we are excluding
    # thus, possible number of HPD intervals
    # this gives us the head of the distribution
    n_intervals = n - interval_idx_inc
    # the width of each possible interval
    # for each possible head and tail value, 
    # what is the difference between them
    interval_width = dat[interval_idx_inc:].reset_index(drop=True) - \
                     dat[:n_intervals].reset_index(drop=True)
    # find the shortest interval
    min_idx = interval_width.values.argmin()
    hpd_interval = (dat.iloc[min_idx], dat.iloc[min_idx+interval_idx_inc])
    dat_hpd = dat[dat.between(hpd_interval[0],
                              hpd_interval[1])]
    dat_mid = dat_hpd.quantile(qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


# reads in beast logs and phydyn trajectory files, calculates HPD
# for parameters in log, returns one dictionary element with 
# all values and one with only values in the HPD. 
# for trajectory items it returns the cumualtive amount (amount at end of time series)
# which fall within the HPD as well as a trajectory withe the median and HPD
# limits for each time step
def process_results(args):
    results = {}
    log = pd.read_csv(args.log, sep="\t", comment="#")
    #plot_trace(log, args)
    log = log[log['Sample']>=max(log['Sample'])*args.burnin]
    for d in log.columns[1:]:   # excludes sample column
        results[d] = log[d]
        limits = hpd(log[d], args.qs)
        results[d+'_HPD'] = log[d][log[d].between(limits[0], 
                                                     limits[2])]
    # reads in PhyDyn trajectories
    # slow but i'm not sure how to make it faster
    traj = pd.read_csv(args.traj, sep='\t')
    traj = traj[traj['Sample']>=max(traj['Sample'])*args.burnin].reset_index(drop=True)
    max_time = max(traj['t'])
    for item in traj.columns[2:]:   # skips time and sample columns
        results = process_traj(results, max_time, traj, item, args)
    # state probabilities (P(new infection from exog))
    state_prob = traj.copy()
    tau = (1-args.ph)*(1/args.ph)/(1/args.P - 1)
    # dictionary with sample as key and R0 as value
    R0_dict = log.set_index('Sample')['seir.R0'].to_dict()
    # sample is key and a is value
    a_dict = log.set_index('Sample')['seir.a'].to_dict()
    state_prob['beta'] = state_prob['Sample'].map(R0_dict)*(366/args.infPeriod)/\
                            ((1-args.ph)+tau*args.ph)
    # accounts for piecewise R0
    state_prob.loc[state_prob['t'] > args.R0ChangeDate, 'beta'] = \
        state_prob['Sample'].map(a_dict) * \
        state_prob.loc[state_prob['t'] > args.R0ChangeDate, 'beta']
    state_prob.loc[state_prob['t'] <= args.sirT0, 'beta'] = 0
    N = max(traj['S'])
    if args.empirical: 
        state_prob['eta'] = 0
        state_prob.loc[state_prob['t'] > args.peakImportDate, 'eta'] = \
            (args.peakImport*np.exp(args.theta*args.importS*(state_prob['t']-args.peakImportDate)))
        state_prob.loc[(args.importT0 < state_prob['t']) & 
                       (args.peakImportDate >= state_prob['t']) , 'eta'] = \
            (np.exp(args.theta*args.importR*(state_prob['t']-args.importT0)))
        state_prob['E_state_prob'] = state_prob['eta']/\
                                     (state_prob['eta'] + state_prob['beta']*\
                                           (state_prob['S']/N)*\
                                           (state_prob['Il']+tau*state_prob['Ih']))
    else: 
        state_prob['E_state_prob'] = args.eta/\
                               (args.eta + state_prob['beta']*\
                                           (state_prob['S']/N)*\
                                           (state_prob['Il']+tau*state_prob['Ih']))
    eta = state_prob[['t', 'eta']]
    eta = pd.DataFrame(state_prob.groupby('t')['eta'].apply(hpd, args.qs).to_list(), 
                 columns=args.qs)
    eta.index = traj['t'][0:len(eta)]
    state_prob = state_prob[['t', 'E_state_prob']]
    state_prob = pd.DataFrame(state_prob.groupby('t')['E_state_prob'].apply(hpd, args.qs).to_list(), 
                 columns=args.qs)
    state_prob.index = traj['t'][0:len(state_prob)]
    results['E_state_prob_HPD'] = state_prob
    results['eta'] = eta
    return(results)


def plot_trace(log, args):
    # remove first 10% as burnin
    # this helps with the axis scales 
    log = log[log['Sample']>=max(log['Sample'])*0.10]
    plot_items = log.columns[1:]
    n_rows = math.ceil(len(plot_items)/2.0)
    fig, axs = plt.subplots(n_rows, 2, 
                            figsize=(6.4*4, 4.8*n_rows),
                            constrained_layout=True)
    for i, item in enumerate(plot_items):
        row = math.floor(i/2)
        col = math.ceil(i/2 - math.floor(i/2))
        axs[row][col].plot(log['Sample'], log[item])
        axs[row][col].set_xlabel('state')
        axs[row][col].set_ylabel(item)
    fig.savefig('figures/'+args.plot_name+'_trace.pdf')
    plt.close(fig)


# Plots trajectory of selected p_h values on the same plot
# TODO make code more concise
def plot_traj_combined(args, all_results, ph_s, eta, israel_reported_data):
    fig, axs = plt.subplots(1, 2, figsize=(6.4*2, 4.8),
                            constrained_layout=True)
    sns.scatterplot(x="numeric_date", y="cumulative cases",
                    data=israel_reported_data,
                    color=plot_style.polar_night[0], ax=axs[1])
    for i, ph in enumerate(ph_s):
        infected_traj = all_results[(ph, eta)]['Il_traj_HPD'] + \
                        all_results[(ph, eta)]['Ih_traj_HPD']
        recovered_traj = all_results[(ph, eta)]['R_traj_HPD']
        axs[0].plot(infected_traj.index,
                    infected_traj[0.5],
                    linewidth=3,
                    color=plot_style.ph_col_dict[ph])
        axs[0].fill_between(infected_traj.index,
                            infected_traj[0.025],
                            infected_traj[0.975],
                            color=plot_style.ph_col_dict[ph],
                            alpha=0.25)
        axs[0].set_ylabel('Infected individuals')
        axs[0].text(0.05, 0.9-i*0.07, '$p_h$ = {0}'.format(ph),
                    color=plot_style.ph_col_dict[ph],
                    size=20, transform=axs[0].transAxes)
        axs[1].plot(recovered_traj.index,
                    recovered_traj[0.5],
                    linewidth=3,
                    color=plot_style.ph_col_dict[ph])
        axs[1].fill_between(recovered_traj.index,
                            recovered_traj[0.025],
                            recovered_traj[0.975],
                            color=plot_style.ph_col_dict[ph],
                            alpha=0.25)
        axs[1].set_ylabel('Cumulative incidence')
    for i, ax in enumerate(axs):
        ax.ticklabel_format(useOffset=False)
        ax.set_yscale('log')
        ax.set_xlabel('Month (2020)')
        ax.set_ylim(1, )
        ax.set_xlim(2020+(31+1)/366-0.035,
                    2020+(31+29+31+31+1)/366+0.035)
        ax.set_xticks((2020+(31+1-0.5)/366, 2020+(31+29+1-0.5)/366,
                       2020+(31+29+31+1-0.5)/366, 2020+(31+29+31+31+1-0.5)/366))
        ax.set_xticklabels(('Feb.', 'March', 'April', 'May'))
        ax.text(-0.15, 1, string.ascii_uppercase[i], color=plot_style.polar_night[0],
                size=24, transform=ax.transAxes)
    fig.savefig('figures/'+args.plot_name+'_traj_combined.pdf')
    plt.close(fig)


# Plots epi parameters across a range of p_h values
def plot_ph_sensitivity(args, data, params, eta, israel_reported_data):
    fig, axs = plt.subplots(1, len(params.keys()), figsize=(6.4*len(params.keys()), 4.8),
                            constrained_layout=True)
    for i, key in enumerate(params.keys()):
        key_name = key.replace('_HPD', '')
        key_name = key_name.replace('_cum', '')
        dat = data[key][data[key]['eta'] == eta]
        # logging the data so that density is calculated on the log scale
        if params[key]['y_scale_log']:
            dat[key_name] = np.log10(dat[key_name])
        sns.violinplot(x='ph', y=key_name, data=dat, 
                       ax=axs[i], saturation=1, cut=0, 
                       palette=plot_style.ph_col_dict)
        axs[i].set_ylabel(params[key]['y_label'])
        axs[i].set_ylim(params[key]['y_lim'])
        axs[i].set_xlabel("$p_h$ (%)")
        if 'y_ticks' in params[key].keys():
            axs[i].set_yticks(params[key]['y_ticks'])
        labels = [item.get_text() for item in axs[i].get_xticklabels()]
        labels = [int(float(item)*100) for item in labels]
        axs[i].set_xticklabels(labels)
        axs[i].text(-0.15, 1, string.ascii_uppercase[i], color=plot_style.polar_night[0],
                size=24, transform=axs[i].transAxes)
        if params[key]['y_scale_log']:
            ticks = [tick for tick in range(0, params[key]['y_lim'][1], 2)]
            axs[i].set_yticks(ticks)
            axs[i].set_yticklabels(
                   ["$\mathregular{{10^{0}}}$".format(x) for x in ticks])
        if params[key]['reported_data']:
            axs[i].axhline(y=np.log10(max(israel_reported_data['cumulative cases'])),
                           linestyle='--', color=plot_style.polar_night[0])
    fig.savefig('figures/'+args.plot_name+'_ph_sensitivity.pdf')
    plt.close(fig)


def plot_ph_eta_sensitivity(args, data, params, israel_reported_data):
    etas = args.etas
    fig, axs = plt.subplots(len(params.keys()), len(etas), 
                            figsize=(6.4*len(etas), 4.8*len(params.keys())), 
                            constrained_layout=True)
    for j, eta in enumerate(etas):
        for i, key in enumerate(params.keys()):
            key_name = key.replace('_HPD', '')
            key_name = key_name.replace('_cum', '')
            dat = data[key][data[key]['eta'] == eta]
            # logging the data so that density is calculated on the log scale
            if params[key]['y_scale_log']:
                dat[key_name] = np.log10(dat[key_name])
            sns.violinplot(x='ph', y=key_name, data=dat, 
                           ax=axs[i][j], saturation=1, cut=0, 
                           palette=plot_style.ph_col_dict)
            axs[i][j].set_ylim(params[key]['y_lim'])
            labels = [item.get_text() for item in axs[i][j].get_xticklabels()]
            labels = [int(float(item)*100) for item in labels]
            axs[i][j].set_xticklabels(labels)
            axs[i][j].set_ylabel('')
            axs[i][j].set_xlabel('')
            axs[i][j].set_title('')
            if 'y_ticks' in params[key].keys():
                axs[i][j].set_yticks(params[key]['y_ticks'])
            axs[i][0].set_ylabel(params[key]['y_label'])
            if params[key]['y_scale_log']:
                ticks = [tick for tick in range(0, params[key]['y_lim'][1], 2)]
                axs[i][j].set_yticks(ticks)
                axs[i][j].set_yticklabels(
                       ["$\mathregular{{10^{0}}}$".format(x) for x in ticks])
            if params[key]['reported_data']:
                axs[i][j].axhline(y=np.log10(max(israel_reported_data['cumulative cases'])),
                               linestyle='--', color=plot_style.polar_night[0])
        axs[-1][j].set_xlabel("$p_h$ (%)")
        axs[0][j].text(-0.15, 1, string.ascii_uppercase[j], color=plot_style.polar_night[0],
                       size=24, transform=axs[0][j].transAxes)
        if args.empirical:
            axs[0][j].set_title("$\Theta$ = {0}".format(eta))
        else:
            axs[0][j].set_title("$\eta$ = {0}".format(eta))
    fig.savefig('figures/{0}_ph_eta_sensitivity.pdf'.format(args.plot_name))
    plt.close(fig)



def plot_state_probs(args, data, ph_s, etas):
    n_col = len(etas)
    fig, axs = plt.subplots(len(ph_s), n_col, 
                            figsize=(6.4*len(etas), 4.8*len(ph_s)), 
                            constrained_layout=True)
    # j is column
    for j, eta in enumerate(etas):
        # i is row
        for i, ph in enumerate(ph_s):
            dat = data[(ph, eta)]['E_state_prob_HPD']
            axs.flat[i*n_col+j].plot(dat.index,
                           dat[0.5],
                           linewidth=3,
                           color=plot_style.ph_col_dict[0.02])
            axs.flat[i*n_col+j].fill_between(dat.index,
                                   dat[0.025],
                                   dat[0.975], 
                                   color=plot_style.ph_col_dict[0.02],
                                   alpha=0.25)
            if args.empirical:
                axs.flat[i*n_col+j].set_title("$\Theta$ = {0}, $p_h$ = {1}%".format(eta, round(ph*100)))
            else:
               axs.flat[i*n_col+j].set_title("$\eta$ = {0}, $p_h$ = {1}%".format(eta, round(ph*100)))
            axs.flat[i*n_col+j].set_xlim(2020+(31+1)/366-0.035,
                    2020+(31+29+31+31+1)/366+0.035)
            axs.flat[i*n_col+j].set_xticks((2020+(31+1-0.5)/366, 2020+(31+29+1-0.5)/366,
                       2020+(31+29+31+1-0.5)/366, 2020+(31+29+31+31+1-0.5)/366))
            axs.flat[i*n_col+j].set_xticklabels(('Feb.', 'March', 'April', 'May'))
            axs.flat[i*n_col+j].set_xlabel('')
            axs.flat[i*n_col+j].set_ylabel('')
            if ph == ph_s[-1]:
                axs.flat[i*n_col+j].set_xlabel('Month (2020)')
            if eta == etas[0]:
                axs.flat[i*n_col+j].set_ylabel('Prob. imported infection')
    fig.savefig('figures/{0}_state_probs.pdf'.format(args.plot_name))


def make_table(args, rows, data, ph_s, etas):
    for eta in etas:
        table = [[rows[row]['label'], rows[row]['prior']] for row in rows.keys()]
        table.insert(0, ['Parameter', 'Prior']) 
        for ph in ph_s:
            table[0].append(f'p_h={ph}')
            for i, key in enumerate(rows.keys()):
                dat = data[(ph, eta)][f'{key}_HPD']
                if rows[key]['sci']:
                    median = "{:.{}e}".format(dat.median(), rows[key]['decimals'])
                    low_limit = "{:.{}e}".format(min(dat), rows[key]['decimals'])
                    high_limit = "{:.{}e}".format(max(dat), rows[key]['decimals'])
                else:
                    median = "{:.{}f}".format(dat.median(), rows[key]['decimals'])
                    low_limit = "{:.{}f}".format(min(dat), rows[key]['decimals'])
                    high_limit = "{:.{}f}".format(max(dat), rows[key]['decimals'])
                table[i+1].append(f'{median} ({low_limit} {high_limit})')
        with open(f'tables/{args.plot_name}_{eta}.csv', 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(table)


def process_reported_data(data_path, max_date):
    reported_data = pd.read_csv(data_path)
    israel_reported_data = reported_data[
        reported_data['countriesAndTerritories'] == 'Israel']
    israel_reported_data['numeric_date'] = pd.to_datetime(
        israel_reported_data['dateRep'], dayfirst=True).apply(numeric_date)
    israel_reported_data.sort_values(by='numeric_date', inplace=True)
    israel_reported_data['cumulative cases'] = israel_reported_data[
        'cases'].cumsum()
    israel_reported_data = israel_reported_data[israel_reported_data[
        'numeric_date'] <= max_date]
    return(israel_reported_data)


def process_results_main(args):
    violin_data = {}
    all_results = {}
    for ph in args.phs:
        for eta in args.etas: 
            print(ph, eta)
            args.ph = ph
            args.eta = eta
            if args.empirical == True:
                args.theta = eta
                args.peakImport = np.exp(args.theta*args.importR*(args.peakImportDate-args.importT0))
            args.log = 'data/{0}/{0}_{1}_{2}/{0}_{1}_{2}.log'.\
                            format(args.run_id, args.ph, args.eta)
            args.traj =  'data/{0}/{0}_{1}_{2}/seir.{0}_{1}_{2}.traj'.\
                            format(args.run_id, args.ph, args.eta)
            args.plot_name = '{0}_{1}_{2}'.format(args.run_id, ph, eta)
            results = process_results(args)
            all_results[(ph, eta)] = results
            for violin in args.violin_data_keys:
                dat = pd.DataFrame(results[violin])
                dat['ph'] = ph
                dat['eta'] = eta
                if violin in violin_data.keys():
                    violin_data[violin] = violin_data[violin].append(dat)
                else:
                    violin_data[violin] = dat
    return(all_results, violin_data)


def plot_reported_data(dat, epi_events):
    fig, axs = plt.subplots(1, 1, 
                            figsize=(6.4, 4.8), 
                            constrained_layout=True)
    axs.plot(dat['numeric_date'], dat['cases'],
        color=plot_style.polar_night[0], lw=3)
    axs.set_xlim((2020.05, 2020.35))
    axs.set_xticks(pd.to_datetime(['02-01-20', '03-01-20', '04-01-20', '05-01-20']).\
                to_series().apply(numeric_date).to_list())
    axs.set_xticklabels(['Feb.', 'March', 'April', 'May'])
    axs.set_xlabel('Month (2020)')
    axs.set_ylabel('Reported cases')
    cols = [plot_style.ph_col_dict[0.02], 
        plot_style.ph_col_dict[0.05], 
        plot_style.ph_col_dict[0.1], 
        plot_style.ph_col_dict[0.2]]
    counter=0
    axs.get_ylim()
    y_pos = 0.95*axs.get_ylim()[1]
    for date, label in epi_events.items():
        axs.axvline(x=date, ls='--', 
            color=cols[counter], 
            alpha=0.75)
        axs.text(date+0.002, y_pos, label, color=cols[counter],
        size=16, verticalalignment='top', 
        bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))
        y_pos = y_pos-0.15*axs.get_ylim()[1]
        counter+=1
    fig.savefig('figures/reported_cases.pdf')
    plt.close('fig')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--log', default='Israel_SARS_CoV-2_constant_R0.log')
    parser.add_argument('--traj', default='seir.Israel_SARS_CoV-2_constant_R0.traj')
    parser.add_argument('--metadata', help='metadata file',
                        default='data/metadata.tsv')
    parser.add_argument('--burnin', default=0.5)
    parser.add_argument('--qs', help='quantiles for summary statistics', default='0.025,0.50,0.975')
    parser.add_argument('--P', help='prop. of infections caused by p_h of population', default=0.8)
    parser.add_argument('--sirT0', help='begining of SIR dynamics', default=(2020+((31+1)/366)))
    parser.add_argument('--R0ChangeDate', help='date on which R0 changes by factor a',
                    default=2020+(31+29+19)/366)
    parser.add_argument('--infPeriod', help='duration of infecitous period (days)',
                    default=5.5)
    parser.add_argument('--run_id', default='')
    parser.add_argument('--empirical', default=False)
    args = parser.parse_args()
    args.qs = tuple([float(item) for item in args.qs.split(',')])
    plot_params = {
                    'seir.R0_HPD': 
                        {'y_label':         r"$ R_0 $",
                         'y_lim':           (0, 4),
                         'y_scale_log':     False,
                         'reported_data':   False},
                    'seir.a_HPD':
                        {'y_label':         r"$ \alpha $",
                         'y_lim':           (0, 2),
                         'y_scale_log':     False,
                         'reported_data':   False,
                         'y_ticks': [0, 0.5, 1.0, 1.5, 2.0]},
                    'R_cum_HPD':     
                        {'y_label':         "Cumulative incidence",
                         'y_lim':           (0, 7),
                         'y_scale_log':     True,
                         'reported_data':   True}}
    table_rows = {  'posterior':    
                        {'label':           'Posterior',
                         'decimals':        0, 
                         'sci':             False,
                         'prior':           ''}, 
                    'TreeHeight':   
                        {'label':           'Tree height',
                         'decimals':        3,
                         'sci':             False,
                         'prior':           ''},
                    'clockRate':    
                        {'label':           'Clock rate',
                         'decimals':        2,
                         'sci':             True,
                         'prior':           f'Uniform(5.0e-4, 2.0e-3)'},
                    'kappa':        
                        {'label':           f'\N{GREEK SMALL LETTER KAPPA}',
                         'decimals':        3,
                         'sci':             False,
                         'prior':           'Lognormal(M=1.0, SD=1.25)'},
                    'gammaShape': 
                        {'label':           f'\N{GREEK SMALL LETTER GAMMA}',
                         'decimals':        2,
                         'sci':             True,
                         'prior':           'Exponential(M=1.0)'},
                    'seir.R0': 
                        {'label':           f'R\N{SUBSCRIPT ZERO}',
                         'decimals':        2,
                         'sci':             False,
                         'prior':           'Lognormal(M=1.5, SD=0.5)'},
                    'seir.a': 
                        {'label':           f'\N{GREEK SMALL LETTER ALPHA}',
                         'decimals':        2,
                         'sci':             False,
                         'prior':           'Uniform(0.0, 2.0)'},
                    'seir.E': 
                        {'label':           f'E\N{SUBSCRIPT ZERO}',
                         'decimals':        2,
                         'sci':             False,
                         'prior':           'Expoenential(M=1.0)'},
                    'seir.exog': 
                        {'label':           f'Y\N{SUBSCRIPT ZERO}',
                         'decimals':        2,
                         'sci':             True,
                         'prior':           'Exponential(M=1.0)'}}
    # Processes results using a constant import rate
    args.phs = [0.02, 0.05, 0.10, 0.20, 0.50, 0.80]
    args.etas = [10, 100, 1000, 2500, 5000]
    args.violin_data_keys = ['seir.R0_HPD', 'seir.a_HPD', 'R_cum_HPD']
    args.run_id = 'june11_CMCMC_Final'
    all_results, violin_data = process_results_main(args)
    with open("{0}_all_results.pkl".format(args.run_id), "wb") as outfile:
            pickle.dump(all_results, outfile)
    with open("{0}_violin_data.pkl".format(args.run_id), "wb") as outfile:
            pickle.dump(violin_data, outfile)
    all_results = pickle.load(open('{0}_all_results.pkl'.format(args.run_id), 'rb'))
    violin_data = pickle.load(open('{0}_violin_data.pkl'.format(args.run_id), 'rb'))
    args.plot_name = 'june11_CMCMC_Final'
    israel_reported_data = process_reported_data('data/ecdc_reported.csv',
                                                 max(all_results[list(
                                                    all_results.keys())[0]]['Il_traj_HPD'].index))
    plot_ph_sensitivity(args, violin_data, plot_params, args.etas[2], israel_reported_data)
    plot_ph_eta_sensitivity(args, violin_data, plot_params, israel_reported_data)
    plot_state_probs(args, all_results, args.phs, args.etas)
    ph_s = [0.02, 0.05, 0.10]
    plot_traj_combined(args, all_results, ph_s, args.etas[2], israel_reported_data)
    make_table(args, table_rows, all_results, args.phs, args.etas)
    # Using empiricial estimates of the import rate
    args.run_id = 'july22_TMRCA_Import_Final'
    args.empirical = True
    args.importT0 = 2020.0511430143656
    args.importR = 56.887337283993205
    args.importS = -51.68986449810366
    args.peakImportDate = 2020.18169399
    all_results = {}
    violin_data = {}
    args.phs = [0.02, 0.05, 0.10, 0.20, 0.50, 0.80]
    args.etas = [0.8, 0.9, 1.0, 1.1, 1.2]
    args.violin_data_keys = ['seir.R0_HPD', 'seir.a_HPD', 'R_cum_HPD']
    all_results, violin_data = process_results_main(args)
    # excluding these parameter values because the mixing is very poor
    for key in violin_data.keys():
        col_name = key.split('_')[0]
        violin_data[key].loc[(violin_data[key]['ph'] == 0.02) & 
            (violin_data[key]['eta'] == 0.8), col_name] = np.nan
        violin_data[key].loc[(violin_data[key]['ph'] == 0.1) & 
            (violin_data[key]['eta'] == 0.9), col_name] = np.nan
    all_results[(0.02, 0.8)]['E_state_prob_HPD'][0.025] = np.nan
    all_results[(0.02, 0.8)]['E_state_prob_HPD'][0.5] = np.nan
    all_results[(0.02, 0.8)]['E_state_prob_HPD'][0.975] = np.nan
    all_results[(0.1, 0.9)]['E_state_prob_HPD'][0.025] = np.nan
    all_results[(0.1, 0.9)]['E_state_prob_HPD'][0.5] = np.nan
    all_results[(0.1, 0.9)]['E_state_prob_HPD'][0.975] = np.nan
    # saves processed data
    with open("{0}_all_results.pkl".format(args.run_id), "wb") as outfile:
            pickle.dump(all_results, outfile)
    with open("{0}_violin_data.pkl".format(args.run_id), "wb") as outfile:
            pickle.dump(violin_data, outfile)
    all_results = pickle.load(open('{0}_all_results.pkl'.format(args.run_id), 'rb'))
    violin_data = pickle.load(open('{0}_violin_data.pkl'.format(args.run_id), 'rb'))
    args.plot_name = 'july22_TMRCA_Import_Final'
    israel_reported_data = process_reported_data('data/mecdc_reported.csv',
                                                 max(all_results[list(
                                                    all_results.keys())[0]]['Il_traj_HPD'].index))
    plot_ph_sensitivity(args, violin_data, plot_params, args.etas[2], israel_reported_data)
    ph_s = [0.02, 0.05, 0.10]
    plot_traj_combined(args, all_results, ph_s, args.etas[2], israel_reported_data)
    plot_ph_eta_sensitivity(args, violin_data, plot_params, israel_reported_data)
    plot_state_probs(args, all_results, args.phs, args.etas)
    make_table(args, table_rows, all_results, args.phs, args.etas)
    args.ph = 0.05
    args.eta = 1.0
    mcc_tree = calc_mcc_tree('data/{0}/{0}_{1}_{2}/{0}_{1}_{2}.trees'.\
                                format(args.run_id, args.ph, args.eta), int(args.burnin*100))
    plot_tree(mcc_tree, args)


if __name__ == "__main__":
    main()


