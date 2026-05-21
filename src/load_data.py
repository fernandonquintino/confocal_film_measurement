import glob

import numpy as np
import pandas as pd


def load_labview(data_path, date, sensor):
    path = f"{data_path}/{date}_{sensor}/"
    csv_files = glob.glob(path + "*.csv")
    li = []

    for filename in csv_files:
        df = pd.read_csv(filename, header=23, encoding='unicode_escape', decimal = ',', sep = '\t').drop(['X_Value', 'Comment', 'DP3'], axis=1)
        li.append(df)

    df = pd.DataFrame()
    for i in range(len(li)):
        df = pd.concat([df, pd.DataFrame(li[i]).mean()], axis=1)
    df = df.T
    df.reset_index(drop=True, inplace=True)

    df.rename(columns={'DP 1': 'Pv', 'DP2': 'dp_exp'}, inplace=True)
    df['sensor'] = sensor

    # df['Qv'] = df['Qv'] - df['Qv'][0]
    df['dp_exp'] = df['dp_exp'] - df['dp_exp'][0]

    df.rename(columns={'T2': 'Tg', 'Qv': 'Qg', 'Pv': 'Pg'}, inplace=True)

    df = df.iloc[1:,:].reset_index(drop=True)

    return df

def load_confocal(data_path, date, header):
    path = f"{data_path}/{date}/"
    csv_files = glob.glob(path + "*.csv")
    li = []

    for filename in csv_files:
        df = pd.read_csv(filename, sep=';', decimal=',', header=header)

        df.rename({'Time stamp - Timestamp': 'time',
                   'Difference 1 / 2 - MeasSignal': 't1',
                   'Difference 2 / 3 - MeasSignal': 't2',
                   'Difference 3 / 4 - MeasSignal': 't3'}, axis=1, inplace=True)

        if 't3' not in df.columns:
            df['t3'] = -10

        df.loc[df['t1']>8, 't1'] = -10
        df.loc[df['t2']>10, 't2'] = -10
        df.loc[df['t3']>8, 't3'] = -10

        df.replace(-10, np.nan, inplace=True)
        df['time'] = df['time'] - df['time'].min()
        df.loc[df['t2']>8.02, 't2'] = 8.02

        li.append(df[['time', 't1', 't2', 't3']])

    return li

def load_confocal2(filename, header):

    df = pd.read_csv(filename, sep=';', decimal=',', header=header)
    df.rename({'Time stamp - Timestamp': 'time',
               'Difference 1 / 2 - MeasSignal': 't1',
               'Difference 2 / 3 - MeasSignal': 't2',
               'Difference 3 / 4 - MeasSignal': 't3'}, axis=1, inplace=True)
    if 't3' not in df.columns:
        df['t3'] = -10
    df.loc[df['t1']>8, 't1'] = -10
    df.loc[df['t2']>10, 't2'] = -10
    df.loc[df['t3']>8, 't3'] = -10

    df.replace(-10, np.nan, inplace=True)
    df['time'] = df['time'] - df['time'].min()
    df.loc[df['t2']>8.02, 't2'] = 8.02

    return df[['time', 't1', 't2', 't3']]
