
import os
# added for later use
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import warnings
warnings.filterwarnings('ignore')
from process_beast_log import *
import argparse


def process_data(log_file):
    f_data = load_data(log_file)
    gr = get_growth_rate(f_data)
    tmrca = tmrca_and_clock_rate(f_data)
    f_data["age(root).calander"] = f_data["age(root)"].apply(numericToCalander)
    f_data["age(root).calander"] = f_data["age(root).calander"].astype("datetime64")
    f_tmrca = f_data.groupby(f_data["age(root).calander"].dt.date).count()
    f_tmrca['datetime'] = pd.to_datetime(f_tmrca.index)
    f_tmrca.set_index('datetime', inplace=True)
    # Ne Tau
    ne_tau_df = ne_tau(f_data)
    # prevelance
    prev, prev_list = skyline_prevalence(f_data)
    # incidence
    incidence = total_incidance(prev_list)
    return f_data, ne_tau_df, prev, incidence

def plot_GR_TMRCA_DT(f_data, out):
    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        curr = f_data.copy()

        curr['doubeling_time'] = curr['exponential.growthRate'].apply(growthRateInYearsToDoubelingTimeInDays)

        sns.kdeplot(curr['exponential.growthRate'], shade=True, alpha=0.2, ax=ax[0], legend=False)
        ax[0].axvline(curr['exponential.growthRate'].quantile(0.025), color='red', alpha=0.5, linestyle='dashed')
        ax[0].axvline(curr['exponential.growthRate'].quantile(0.975), color='red', alpha=0.5, linestyle='dashed')
        ax[0].axvline(curr['exponential.growthRate'].median(), color='red', alpha=0.5)
        ax[0].set_title("Growth rate")

        sns.kdeplot(curr['doubeling_time'], shade=True, alpha=0.2, ax=ax[1], legend=False)
        ax[1].axvline(curr['doubeling_time'].quantile(0.025), color='red', alpha=0.5, linestyle='dashed')
        ax[1].axvline(curr['doubeling_time'].quantile(0.975), color='red', alpha=0.5, linestyle='dashed')
        ax[1].axvline(curr['doubeling_time'].median(), color='red', alpha=0.5)
        ax[1].set_title("Doubeling time")

        sns.kdeplot(curr['age(root)'], shade=True, alpha=0.2, ax=ax[2], legend=False)
        ax[2].axvline(curr['age(root)'].quantile(0.025), color='red', alpha=0.5, linestyle='dashed')
        ax[2].axvline(curr['age(root)'].quantile(0.975), color='red', alpha=0.5, linestyle='dashed')
        ax[2].axvline(curr['age(root)'].median(), color='red', alpha=0.5)
        ax[2].set_title("tMRCA")

        plt.tight_layout()

        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.clf()
        plt.cla()

def plot_neTau(ne_tau_df, out):
    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        fig, ax = plt.subplots(figsize=(12, 4))
        curr = ne_tau_df
        #curr = curr[curr['time'] >= 2020.12]  # remove dates befor start of february
        sns.lineplot(x='time_calander', y='Ne_tau_median', data=curr, ax=ax, color='red', alpha=0.6)
        ax.fill_between(curr["time_calander"], curr["Ne_tau_upper"],
                        curr["Ne_tau_lower"], alpha=0.2, color='#E38074')
        plt.yscale('log')
        plt.xticks(rotation=90)
        plt.ylabel(r'$N_e\tau$', fontsize=18)
        plt.xlabel('Time')
        sns.despine()
        plt.tight_layout()

        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.clf()
        plt.cla()


def plot_point_prevalence(prev, out):
    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        curr = prev
        #curr = curr[curr['time'] >= 2020.12]  # remove dates befor start of february
        sns.lineplot(x='time_calander', y='prev_median', data=curr, color='red', alpha=0.6)
        plt.fill_between(curr["time_calander"], curr["prev_upper"],
                         curr["prev_lower"], alpha=0.3, color='#E38074')
        plt.yscale('log')
        plt.xticks(rotation=90)
        plt.ylabel('Point prevalence', fontsize=18)
        plt.xlabel('Time')
        sns.despine()
        plt.tight_layout()

        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.clf()
        plt.cla()

def plot_total_incidence(incidence, out):
    with sns.plotting_context(rc={"font.size": 14, "axes.titlesize": 18, "axes.labelsize": 18,
                                  "xtick.labelsize": 14, "ytick.labelsize": 14, 'y.labelsize': 16}):
        curr = incidence
        # curr = curr[curr['time'] >= 2020.12] # remove dates befor start of february
        sns.lineplot(x='time_calander', y='incidence_median', data=curr, color='green', alpha=0.6)
        plt.fill_between(curr["time_calander"], curr["incidence_upper"],
                         curr["incidence_lower"], alpha=0.3, color='#387827')
        plt.yscale('log')
        plt.xticks(rotation=90)
        plt.ylabel('Total incidence', fontsize=18)
        plt.xlabel('Time')
        sns.despine()
        plt.tight_layout()

        if out != None:
            plt.savefig(out, format='png', dpi=350, layout='tight')
        plt.clf()
        plt.cla()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="log file full path", required=True)
    parser.add_argument("-o", "--out", type=str, help="output directory", required=True)
    parser.add_argument("-m", "--mode", type=int, help="the mode of running: 1 regular, 2 partition 3 taxon", default=1)
    args = parser.parse_args()

    log_file = args.file
    mode = args.mode
    out = args.out
    if mode == 1:
        log_files = [log_file]
    elif mode == 2:
        log_files = split_log_with_partitions(log_file)
    elif mode == 3:
        log_files = [split_log_with_taxa(log_file)]
    else:
        raise Exception("running mode is not valid, should be 1 (regular) ,2 (partition) or 3 (taxon)")

    for f in tqdm(log_files):
        print(f"Processing {len(log_files)} log files...")
        f_data, ne_tau_df, prev, incidence = process_data(f)
        base = os.path.basename(f).split('.log')[0]
        plot_GR_TMRCA_DT(f_data, out=os.path.join(out, f"Params_{base}.png"))
        plot_neTau(ne_tau_df, out=os.path.join(out, f"NeTau_{base}.png"))
        plot_point_prevalence(prev, out=os.path.join(out, f"Prevalence_{base}.png"))
        plot_total_incidence(incidence, out=os.path.join(out, f"Incidence_{base}.png"))












