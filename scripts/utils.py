import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from datetime import datetime as dt
import time
import json

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def process_fasta_to_phydyn(sequences, metadata, output):
    records = []
    for rec in SeqIO.parse(sequences, 'fasta'):
        date = pd.to_datetime(metadata[metadata['strain'] == rec.id]['date'].values[0])
        if 'Israel' in rec.id:
            rec.id = rec.id.replace('/', '-').replace('_','-') + \
                     '_' + str(toYearFraction(date)) + '_Il'
        else:
            rec.id = rec.id.replace('/', '-').replace('_','-') + \
                     '_' + str(toYearFraction(date)) + '_exog'
        rec.name = ''
        rec.description = ''
        records.append(rec)
    cnt = SeqIO.write(records, output, 'fasta')
    return cnt


def flatten_json(nested_json):
    """
        Flatten json object with nested keys into a single level.
        Args:
            nested_json: A nested json object.
        Returns:
            The flattened json object if successful, None otherwise.
    """
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(nested_json)
    return out

def get_movement_data(ncov):
    """
    TODO adjust code - need its own script?
    :param ncov:
    :return:
    """
    flattened = flatten_json(ncov['tree'])
    isr_leaves = {k:v for k,v in flattened.items() if 'Israel/' in str(v) and '_name' in k}
    df = pd.DataFrame(isr_leaves.items(), columns=['id', 'tip'])
    df['group1'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-3]))
    df['group2'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-5]))
    df['group3'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-7]))
    df['group4'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-9]))
    df['group5'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-11]))
    df['group6'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-13]))
    # df['group7'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-15]))
    # df['group8'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-17]))
    try:
        df['country_anc1'] = df['group1'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        df['country_anc2'] = df['group2'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        df['country_anc3'] = df['group3'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        df['country_anc4'] = df['group4'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        df['country_anc5'] = df['group5'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        df['country_anc6'] = df['group6'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        # df['country_anc7'] = df['group7'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
        # df['country_anc8'] = df['group8'].apply(lambda val: flattened[val + "_node_attrs_country_value"])
    except:
        print("ancestry retrieval error, num columns: {}".format(len(df.columns)))

    df['id2'] = df['id'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    df['country_self'] = df['id2'].apply(lambda val: flattened[val + "_node_attrs_country_value"])

    # get first ancestral country change
    df['ancestor_first_change'] = None
    df['country_first_change'] = None

    anc_count = (len(df.columns) - 4) / 2
    print('anc count: {}'.format(anc_count))

    for i in reversed(range(1,anc_count + 1)):
        df['ancestor_first_change'] = np.where(df['country_anc{}'.format(i) ] != df['country_self'], i, df['ancestor_first_change'])
        df['country_first_change'] = np.where(df['country_anc{}'.format(i) ] != df['country_self'], df['country_anc{}'.format(i) ], df['country_first_change'])

    df['hierarcy_first_change'] = None
    for index, row in df.iterrows():
        if row['ancestor_first_change'] is not None:
            group_first_change = 'group'+ str(row['ancestor_first_change'])
            row['hierarcy_first_change'] = row[group_first_change]

    # all_movements_per_tip
    all_movements_per_tip = df[['tip', 'country_self'
                                    ,'country_anc1'
                                    , 'country_first_change'
                                    , 'hierarcy_first_change'
                                   ]]

    all_movements_per_tip = all_movements_per_tip.sort_values(by=['country_self', 'country_anc1', 'country_first_change'])

    # in-israel movements
    internal_only_movements_per_tip = all_movements_per_tip[
        (all_movements_per_tip['country_first_change'].isin(
            ['North District', 'South District', 'Tel Aviv District', 'Jerusalem District', 'South Coast District']))
    ]

    # movements per region
    # jrsm_movements = all_movements_per_tip[
    #     (all_movements_per_tip['country_first_change'].isin(
    #         ['Jerusalem']))
    # ]
    # jrsm_movements.sort_values(by=['country_self', 'country_anc1', 'country_first_change'], inplace=True)

    return all_movements_per_tip, internal_only_movements_per_tip


def get_exposure_statistics(movement_data):

    # verification - get #unique movements according to position on tree
    unique_verification = movement_data.groupby(['hierarcy_first_change']).size().reset_index(name='size')
    unique_verification_count = unique_verification.count()
    print(unique_verification_count)

    # external - summarize unique movements by source country
    unique_source = movement_data.groupby(['country_first_change', 'hierarcy_first_change']).size().reset_index(name='size')
    unique_source_count = unique_source.groupby('country_first_change').size().reset_index(name='size')
    unique_source_count.sort_values(by=['size'], ascending=False, inplace=True)

    # internal - summarize unique movements by source+target country
    unique_source_target = movement_data.groupby(['country_first_change', 'country_self', 'hierarcy_first_change']).size().reset_index(name='size')
    unique_source_target_count = unique_source_target.groupby(['country_first_change', 'country_self']).size().reset_index(name='size')
    unique_source_target_count.sort_values(by=['size'], ascending=False, inplace=True)


    return unique_source_count, unique_source_target_count


def main_internal():
    ncov_json_path = '/Users/omer/Documents/Lab/COVID/auspice/ncov_districts.json'
    with open(ncov_json_path) as json_file:
        ncov = json.load(json_file)

    all_movement_data, internal_only_movement_data = get_movement_data(ncov)
    # export to files
    # all_movement_data.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/all_movement_data.csv', index=False)
    internal_only_movement_data.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/internal_only_movement_data.csv', index=False)

    exposure_summary_external, exposure_summary_internal = get_exposure_statistics(all_movement_data)
    # export to files
    # exposure_summary_external.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/exposure_summary_external.csv', index=False)
    exposure_summary_internal.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/exposure_summary_internal.csv', index=False)


def main_external_batch(dir_name):

    files = glob.glob('/sternadi/home/volume3/COVID19/data/bootsrap/{}/*/results/ncov_with_accessions_and_travel_branches.json'.format(dir_name))

    df = None
    for f in files:
        print(f)
        run_id = (f.split('/')[8]).split('_')[1]
        with open(f) as json_file:
            ncov = json.load(json_file)

        all_movement_data, internal_only_movement_data = get_movement_data(ncov)
        # export to files
        # all_movement_data.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/all_movement_data.csv', index=False)
        # internal_only_movement_data.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/internal_only_movement_data.csv', index=False)

        exposure_summary_external, exposure_summary_internal = get_exposure_statistics(all_movement_data)
        # export to files
        # exposure_summary_external.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/exposure_summary_external.csv', index=False)
        # exposure_summary_internal.to_csv('/Users/omer/Documents/Lab/COVID/Analysis/exposure_summary_internal.csv', index=False)

        exposure_summary_external['size'] = exposure_summary_external['size'] / exposure_summary_external['size'].sum()

        if df is None:
            df = exposure_summary_external
        else:
            df = df.merge(exposure_summary_external, how='outer', on='country_first_change', sort= True, suffixes=('', '_'+str(run_id)))

    df.to_csv('/sternadi/home/volume3/omer/repos/COVID19/results/exposure_{}.csv'.format(dir_name), index=False)


if __name__ == '__main__':
    main_external_batch('world')


