import pandas as pd
import PyAstronomy.pyasl
from datetime import date
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import random
from tqdm import tqdm
import re

#Auxilary functions
def growthRateInYearsToDoubelingTimeInDays(g):
    return (365 * (np.log(2) / g))


def numericToCalander(d):
    return PyAstronomy.pyasl.decimalYearGregorianDate(d)


def create_skylineNeTau_plotData(tmrca, popSize, growthRate, starDate=2019.9, endDate=2020.25, interval=0.005):
    times = np.arange(starDate, endDate, interval)
    skylineNeTau = {t:np.exp(growthRate*(t-tmrca))/popSize for t in times}
    return pd.DataFrame(skylineNeTau.items(), columns=["time", "Ne_tau"])


def getSigmaSquared(meanR0, dispersion):
    return((1/dispersion) + (1/meanR0) + 1)


def create_skylinePrev_plotData(tmrca, popSize, growthRate, starDate=2019.9, endDate=2020.25, interval=0.005):
    meanR0 = random.uniform(1.8, 2.8)
    k = random.uniform(0.15, 0.3)
    """
    NOTE - in TB he samples random numbers for meanR0 and K
    here I (Talia) chose specific parameters and also chose and ran them in TBs notebook 
    so I could compare the results
    When we are done editing we can remove the next 2 lines
    """
    meanR0 = 2.4
    k = 0.16
    tau = 0.020
    sigmaSquared = getSigmaSquared(meanR0, k)
    times = np.arange(starDate, endDate, interval)
    skylinePrev = {t: (meanR0 * (1 + 1/k)) * (np.exp(growthRate * (t - tmrca)) / (popSize * tau)) for t
                   in times}
    return pd.DataFrame(skylinePrev.items(), columns=["time", "Prevalence"])


def get_skylineNeTau_atTime(skylineNeTauList, time):
    deltaT = 0.004
    return (skylineNeTauList.loc[(skylineNeTauList.time >= time-deltaT) & (skylineNeTauList.time < time+deltaT)].drop("time", axis=1))


def quantileSkyline(skylineNeTauList, quantile, starDate=2019.9, endDate=2020.25, interval=0.005, colname=""):
    df = pd.DataFrame()
    times = np.arange(starDate, endDate, interval)
    for t in times:
        df = df.append({"time":t, colname:np.quantile(get_skylineNeTau_atTime(skylineNeTauList, t), quantile)}, ignore_index=True)
    return(df)


def load_data(input_file):
    try:
        df = pd.read_csv(input_file, sep="\t") # sometimes we change the file - no 3 rows header (partitions or taxon sets)
    except:
        df = pd.read_csv(input_file, skiprows=3, sep="\t")
    # parameters from TBs notebook - we can change it as we want
    burnin = 100000000 * 0.4
    skip = 2
    data = df.loc[df.state >= burnin]
    data = data.iloc[::skip, :]
    return data

def get_growth_rate(data):
    # calculate the posterior from the MCMC
    growthRate = data["exponential.growthRate"].median()
    growthRateLowerQuantile = data["exponential.growthRate"].quantile(0.025)
    growthRateUpperQuantile = data["exponential.growthRate"].quantile(0.975)
    doubelingTime = growthRateInYearsToDoubelingTimeInDays(growthRate)
    doubelingTimeLowerQuantile = growthRateInYearsToDoubelingTimeInDays(growthRateLowerQuantile)
    doubelingTimeUpperQuantile = growthRateInYearsToDoubelingTimeInDays(growthRateUpperQuantile)
    print(
        f"Growth rate: {np.round(growthRate, 2)} [95% CI {np.round(growthRateLowerQuantile, 2)}-{np.round(growthRateUpperQuantile, 2)}]")
    print(
        f"Doubling time (days): {np.round(doubelingTime, 2)} [95% CI {np.round(doubelingTimeUpperQuantile, 2)}-{np.round(doubelingTimeLowerQuantile, 2)}]")

    vals = (growthRate, growthRateLowerQuantile, growthRateUpperQuantile, doubelingTime,
            doubelingTimeLowerQuantile, doubelingTimeUpperQuantile)
    return vals


def tmrca_and_clock_rate(data):
    rootAge = numericToCalander(data["age(root)"].median())
    rootAgeLowerQuantile = numericToCalander(data["age(root)"].quantile(0.025))
    rootAgeUpperQuantile = numericToCalander(data["age(root)"].quantile(0.975))
    clockRate = np.round(data["clock.rate"].median(), 6)
    clockRateLowerQuantile = np.round(data["clock.rate"].quantile(0.025), 6)
    clockRateUpperQuantile = np.round(data["clock.rate"].quantile(0.975), 6)

    print(f"Root age: {rootAge} [95% CI {rootAgeLowerQuantile} to {rootAgeUpperQuantile}]")

    vals = (clockRate, clockRateLowerQuantile, clockRateUpperQuantile, rootAge,
            rootAgeLowerQuantile, rootAgeUpperQuantile)
    return vals


def ne_tau(data):
    skylineNeTauList = pd.DataFrame()
    for index, row in tqdm(data.iterrows()):
        skylineNeTauList = skylineNeTauList.append(
            create_skylineNeTau_plotData(row["age(root)"], row["exponential.popSize"], row["exponential.growthRate"]))
    medianNeTau = quantileSkyline(skylineNeTauList, 0.5, colname="Ne_tau_median")
    lowerNeTau = quantileSkyline(skylineNeTauList, 0.025, colname="Ne_tau_lower")
    upperNeTau = quantileSkyline(skylineNeTauList, 0.975, colname="Ne_tau_upper")
    NeTauForPlot = medianNeTau.merge(lowerNeTau, on="time").merge(upperNeTau, on="time")
    NeTauForPlot["time_calander"] = NeTauForPlot["time"].apply(numericToCalander)
    return NeTauForPlot


def skyline_prevalence(data):
    skylinePrevList = pd.DataFrame()
    for index, row in tqdm(data.iterrows()):
        skylinePrevList = skylinePrevList.append(
            create_skylinePrev_plotData(row["age(root)"], row["exponential.popSize"], row["exponential.growthRate"]))
    medianPrev = quantileSkyline(skylinePrevList, 0.5, colname="prev_median")
    lowerPrev = quantileSkyline(skylinePrevList, 0.025, colname="prev_lower")
    upperPrev = quantileSkyline(skylinePrevList, 0.975, colname="prev_upper")
    prevForPlot = medianPrev.merge(lowerPrev, on="time").merge(upperPrev, on="time")
    prevForPlot["time_calander"] = prevForPlot["time"].apply(numericToCalander)
    return prevForPlot, skylinePrevList


def total_incidance(prevalence):
    times = np.arange(prevalence["time"].min(), prevalence["time"].max() + 0.02, 0.02)

    incidence = pd.DataFrame()
    count_median = 0
    count_lower = 0
    count_upper = 0
    for t in times[:-1]:
        count_median += np.quantile(get_skylineNeTau_atTime(prevalence, t), 0.5)
        count_lower += np.quantile(get_skylineNeTau_atTime(prevalence, t), 0.025)
        count_upper += np.quantile(get_skylineNeTau_atTime(prevalence, t), 0.975)

        incidence = incidence.append({"time": t, "incidence_median": count_median, "incidence_lower": count_lower,
                                      "incidence_upper": count_upper}, ignore_index=True)
    incidence["time_calander"] = incidence["time"].apply(numericToCalander)
    return incidence

# notice- hacky code
def split_log_with_partitions(log_file):

    df = pd.read_csv(log_file, skiprows=3, sep="\t")
    partitions = [c.split('.')[-2] for c in df.columns if 'age(root' in c and 'partition' in c]
    all_log_files = [] # for later use
    for p in partitions:
        df = pd.read_csv(log_file, skiprows=3, sep="\t")
        alias = p
        filtered_columns = [col for col in df.columns if p in col or 'partition' not in col]
        df = df.loc[:, df.columns.isin(filtered_columns)]
        for col in df.columns:
            if p in col:
                if 'treeLikelihood' in p:
                    df.rename(columns={col: 'treeLikelihood'}, inplace=True)    # cannot make every line pretty
                df.rename(columns={col: col.split(f"{p}.")[-1]}, inplace=True)

        new_log_file = log_file.split(".")[0] + '.' + alias + '.log'
        df.to_csv(new_log_file, sep='\t', index=False)
        all_log_files.append(new_log_file)
    return all_log_files


# notice- hacky code
def split_log_with_taxa(log_file):

    df = pd.read_csv(log_file, skiprows=3, sep="\t")
    taxon_name = [re.search('age\((.[a-zA-Z\s]*)', c).group(1) for c in df.columns if 'root' not in c and 'age' in c][0]

    df = df.drop(columns= 'age(root)')
    df.rename(columns={f'age({taxon_name})': 'age(root)'}, inplace=True)

    new_log_file = log_file.split(".log")[0] + '.taxon_' + taxon_name + '.log'
    df.to_csv(new_log_file, sep='\t', index=False)
    return new_log_file

