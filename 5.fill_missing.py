import os
import numpy as np
import pandas as pd

from sklearn.impute import KNNImputer
from sklearn.preprocessing import OrdinalEncoder

DATA_PATH = './variants/SNP_P001.vcf'
rs_max_missing = 30
sample_max_missing = 40

def fill(arr):
    imputer = KNNImputer(n_neighbors=10)
    d = imputer.fit_transform(arr)
    return d


def fill_numerical(Data):
    filled = fill(Data.copy())
    return pd.DataFrame(filled, columns=Data.columns)

def fill_categorical(Cat_data):
    encoder = OrdinalEncoder()

    encoded_data = encoder.fit_transform(Cat_data.values)

    encoded_data[encoded_data == -1] = np.nan
    data = fill(encoded_data)

    encoded_data = np.round(data)

    # Inverse transform the encoded data to get the original data
    df_decoded = pd.DataFrame(encoder.inverse_transform(encoded_data), columns=Cat_data.columns, index=Cat_data.index)

    return df_decoded


def read_vcf(data_path, raw_vcf=True):
    df = pd.read_csv(data_path, delimiter=' ')
    id_col = df.iloc[:,0] + '_' + df.iloc[:, 1].astype(str)   
    cols = df.columns.to_list()[9:]
    print(cols)

    df = df[['[3]ID',*cols]]
    df.ID = id_col

    def trans_col(col):
        try:
            col = col.split('/')[-1].split('_')[0]
        except:
            print(col)
        return col


    cols = [trans_col(col) for col in df.columns]
    df.columns = cols
    df = df.T.rename(columns=df.T.iloc[0]).drop(df.T.index[0]) 
    df = df.replace(['./.'], [np.nan]).apply(lambda x: x.str.replace('/', ''))

    df.to_csv('./variants/original.csv')

    return df


df = read_vcf(DATA_PATH)

missing_rate_samples = ((df.isna().mean(axis=1) * 100).round(2)).sort_values(ascending=False)
missing_rate_snps = ((df.isna().mean(axis=0) * 100).round(2)).sort_values(ascending=False)

cols = missing_rate_snps[missing_rate_snps<rs_max_missing].index
df = df.loc[:, cols[::-1]]

rows = list(missing_rate_samples[missing_rate_samples<sample_max_missing].index)
if 'Undetermined' in rows:
    rows.remove('Undetermined')

print('Following samples are removed:\n',*df.index[~df.index.isin(rows)], sep='\n\t')
df = df.loc[rows[::-1], :]

filled = fill_categorical(df)



def check_alleles(x):
    l = len(set(list(''.join(x))))
    if l >= 3:
        return x.name


war_rs = filled.apply(check_alleles).unique()[1:]
filled.drop(columns=war_rs).to_csv('./variants/filled.csv')
print('\nWARNING!\nFollowing columns have more than 3 alleles and they have been removed from final version of data:',*war_rs, sep='\n\t')

print('Total number of samples:', filled.shape[0])
print('Total number of SNPs:', filled.shape[1])
