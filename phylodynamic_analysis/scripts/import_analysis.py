import pandas as pd
from calendar import isleap
import numpy as np
from scipy.stats import poisson
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sars_cov2_plotting_style as plot_style
import argparse


def numeric_date(dt):
    days_in_year = 366 if isleap(dt.year) else 365
    res = dt.year + (dt.timetuple().tm_yday-0.5) / days_in_year
    return res


def process_data(import_date_file, date_min, date_max):
    import_dates = pd.read_csv(import_date_file)
    # remove rogue January and May date
    import_dates['date_first_change'] = \
        pd.to_datetime(import_dates['date_first_change']).apply(numeric_date).round(8)
    import_dates.columns = ['date', 'importations']
    date_range = np.round(np.arange(date_min, date_max, 1/366), 8)
    # dates not in the data
    missing_dates = [item for item in date_range if item not in list(import_dates['date'])]
    import_dates = import_dates.append(pd.DataFrame.from_dict({import_dates.columns[0]: missing_dates,
    import_dates.columns[1]: [0 for item in missing_dates]}))
    import_dates = import_dates.set_index('date').sort_index()
    peak = import_dates[import_dates['importations'] == max(import_dates['importations'])]
    return(import_dates, peak)


def calc_llh(params, chunked_data, up_times, down_times, chunked_peak_date, tstep):
    r, s, t0 = params
    model_up = np.exp(r*(up_times - t0))
    max_importations = max(model_up)
    model_down = max_importations*np.exp(s*(down_times - chunked_peak_date))
    model = np.append(model_up, model_down)
    model_cum = model.cumsum()*tstep
    chunked_data['model'] = chunked_data.apply(lambda k:
        model_cum[int(k['end_index'])] - model_cum[int(k['start_index'])], axis=1)
    chunked_data['lh'] = poisson.pmf(chunked_data['importations'],
                                             chunked_data['model'])
    chunked_data['llh'] = np.log(chunked_data['lh'])
    return(-sum(chunked_data['llh']))


def plot_importations(chunked_data, chunked_peak_date, up_times, down_times,
                      r, s, import_t0):
    model_up = np.exp(r*(up_times - import_t0))
    model_down = max(model_up)*np.exp(s*(down_times - chunked_peak_date))
    model = np.append(model_up, model_down)
    times = np.append(up_times, down_times)
    fig, axs = plt.subplots(1, 1, figsize=(6.4, 4.8), constrained_layout=True)
    axs.bar((chunked_data['start_time']+chunked_data['end_time'])/2,
            chunked_data['importations']/3,
            3/366, color=plot_style.country_col_dict['Africa'], alpha=0.4)
    axs.plot(times, model*1/366, linewidth=3,
             color=plot_style.country_col_dict['Israel'])
    axs.ticklabel_format(useOffset=False)
    axs.set_xticks((2020+1/366, 2020+(31+1)/366, 2020+(31+29+1)/366,
                   2020+(31+29+31+1)/366, 2020+(31+29+31+31+1)/366))
    axs.set_xticklabels(('Jan.', 'Feb.', 'March', 'April', 'May'))
    axs.set_ylabel('Importations per day')
    axs.set_xlabel('Month (2020)')
    axs.text(0.05, 0.9, 'Model fit',
             color=plot_style.country_col_dict['Israel'],
             size=20, transform=axs.transAxes)
    axs.text(0.05, 0.9-0.07, 'Observed importations',
             color=plot_style.country_col_dict['Africa'],
             size=20, transform=axs.transAxes)
    fig.savefig('figures/importation_dates.pdf')
    plt.close(fig)


def calc_importation_parameters(import_date_file):
    date_min = 2020+0.5/366
    date_max = 2020+(31+29+31+31)/366
    time_step = 3
    continuous_times, tstep = \
        np.linspace(date_min, date_max, 100000, retstep=True)
    chunked_times = np.round(np.arange(date_min, date_max, time_step/366), 8)
    counted_dates, peak = process_data(import_date_file, date_min, date_max)
    chunked_data = counted_dates['importations'].\
        groupby(np.arange(len(counted_dates))//3).sum().to_frame()
    # Counts DO include imoprtations which occur AT start time
    # Counts DO NOT INCLUDE importations which occured AT end time
    chunked_data['start_time'] = np.round(chunked_times, 8)
    chunked_data['end_time'] = np.round(chunked_times+time_step/366, 8)
    chunked_data['start_index'] = \
        chunked_data['start_time'].apply(lambda k:
                                         len(continuous_times[
                                             continuous_times <= k]))
    chunked_data['end_index'] = \
        chunked_data['end_time'].apply(lambda k:
                                       len(continuous_times) -
                                       len(continuous_times[
                                           continuous_times >= k])-1)
    chunked_peak_date = chunked_data[chunked_data['importations'] == max(chunked_data['importations'])]['end_time'].iloc[0]
    up_times = continuous_times[continuous_times <
                                chunked_peak_date]
    down_times = continuous_times[continuous_times >=
                                  chunked_peak_date]
    llh_args = (chunked_data, up_times, down_times, chunked_peak_date, tstep)
    chunked_data
    minimized = minimize(calc_llh, (50, -50, 2020.1), llh_args,
                         method='Nelder-Mead', tol=1e-10,
                         bounds=((1,), (), (2020, 2021)))
    r, s, import_t0 = minimized['x']
    max_up = np.exp(r*(up_times - import_t0))[-1]
    plot_importations(chunked_data, chunked_peak_date, up_times, down_times,
                      r, s, import_t0)
    return(chunked_peak_date, r, max_up, s, import_t0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--importDateFile',
                        default='data/mike/dates_tmrca.csv')
    args = parser.parse_args()
    chunked_peak_date, r, s, import_t0 = calc_importation_parameters(args.importDateFile)
    return(chunked_peak_date, r, s, import_t0)


if __name__ == "__main__":
    main()


