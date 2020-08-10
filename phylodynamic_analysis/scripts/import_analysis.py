import pandas as pd
from calendar import isleap
import numpy as np
from scipy.stats import poisson
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sars_cov2_plotting_style as plot_style
import argparse
import string


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


def plot_importations(time_step, chunked_data, chunked_peak_date, up_times, down_times,
                      r, s, import_t0, multiplier):
    cols = ['#384d67', '#4b6789', '#5e81ac', '#7e9abc', '#8ea6c4']
    up_times = up_times[np.where(up_times > import_t0)]
    model_up = np.exp(r*(up_times - import_t0))
    model_down = max(model_up)*np.exp(s*(down_times - chunked_peak_date))
    model = np.append(model_up, model_down)
    times = np.append(up_times, down_times)
    fig, axs = plt.subplots(1, 2, figsize=(6.4*2, 4.8), constrained_layout=True)
    axs[0].bar((chunked_data['start_time']+chunked_data['end_time'])/2,
            chunked_data['importations']/time_step,
            time_step/366, color=plot_style.country_col_dict['Africa'], alpha=0.4)
    axs[0].plot(times, model*1/366, linewidth=3,
             color=plot_style.country_col_dict['Israel'])
    axs[0].text(-0.15, 1, string.ascii_uppercase[0], color=plot_style.polar_night[0],
                size=24, transform=axs[0].transAxes, fontweight='bold')
    axs[0].ticklabel_format(useOffset=False)
    axs[0].set_xticks((2020+1/366, 2020+(31+1)/366, 2020+(31+29+1)/366,
                   2020+(31+29+31+1)/366, 2020+(31+29+31+31+1)/366))
    axs[0].set_xticklabels(('Jan.', 'Feb.', 'March', 'April', 'May'))
    axs[0].set_ylabel('Importations per day')
    axs[0].set_xlabel('Month (2020)')
    axs[0].text(0.05, 0.9, 'Model fit',
             color=plot_style.country_col_dict['Israel'],
             size=20, transform=axs[0].transAxes)
    axs[0].text(0.05, 0.9-0.07, 'Observed importations',
             color=plot_style.country_col_dict['Africa'],
             size=20, transform=axs[0].transAxes)
    for index, i in enumerate(sorted(multiplier, reverse=True)):
        model_up = np.exp(i*r*(up_times - import_t0))
        model_down = max(model_up)*np.exp(i*s*(down_times - chunked_peak_date))
        model = np.append(model_up, model_down)
        axs[1].plot(times, model*1/366, linewidth=3, color=cols[index])
        axs[1].text(0.05, 0.9-index*0.07, '$\\theta$ = '+str(i),
                 color=cols[index],
                 size=20, transform=axs[1].transAxes)
    axs[1].text(-0.15, 1, string.ascii_uppercase[1], color=plot_style.polar_night[0],
                size=24, transform=axs[1].transAxes, fontweight='bold')
    axs[1].set_ylabel('')
    axs[1].set_ylim(bottom=0)
    axs[1].ticklabel_format(useOffset=False)
    axs[1].set_xticks((2020+1/366, 2020+(31+1)/366, 2020+(31+29+1)/366,
                   2020+(31+29+31+1)/366, 2020+(31+29+31+31+1)/366))
    axs[1].set_xticklabels(('Jan.', 'Feb.', 'March', 'April', 'May'))
    axs[1].set_ylabel('Importations per day')
    axs[1].set_xlabel('Month (2020)')
    fig.savefig(f'figures/importation_dates_{time_step}.pdf')
    plt.close(fig)


def calc_importation_parameters(import_date_file, time_step, multiplier):
    date_min = 2020+0.5/366
    date_max = 2020+(31+29+31+31)/366
    continuous_times, tstep = \
        np.linspace(date_min, date_max, 100000, retstep=True)
    chunked_times = np.round(np.arange(date_min, date_max, time_step/366), 8)
    counted_dates, peak = process_data(import_date_file, date_min, date_max)
    chunked_data = counted_dates['importations'].\
        groupby(np.arange(len(counted_dates))//time_step).sum().to_frame()
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
    minimized = minimize(calc_llh, (50, -50, 2020.1), llh_args,
                         method='Nelder-Mead', tol=1e-10,
                         bounds=((1,), (), (2020, 2021)))
    r, s, import_t0 = minimized['x']
    plot_importations(time_step, chunked_data, chunked_peak_date, up_times, down_times,
                      r, s, import_t0, multiplier)
    import_rates = {}
    for m in multiplier:
        max_up = max(np.exp(m*r*(up_times-import_t0)))
        import_rates[m] = f'if (t > {chunked_peak_date}) then ({max_up}*exp({m}*{s}*(t-{chunked_peak_date}))) else if (t > {import_t0}) then (exp({m}*{r}*(t-{import_t0}))) else 0.0'
    import_rates = pd.DataFrame({'theta': list(import_rates.keys()), 'eta': list(import_rates.values())})
    import_rates.to_csv('data/import_rates.csv',index=None)
    return(chunked_peak_date, r, max_up, s, import_t0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--importDateFile',
        default='data/dates_tmrca.csv')
    parser.add_argument('--tStep', type=float,
        default=3.0)
    parser.add_argument('--multiplier', type=float, nargs='+',
        help='multiplier values', default=[0.8, 0.9, 1, 1.1, 1.2])
    args = parser.parse_args()
    chunked_peak_date, r, max_up, s, import_t0 = calc_importation_parameters(args.importDateFile, args.tStep, args.multiplier)
    
    
if __name__ == "__main__":
    main()


