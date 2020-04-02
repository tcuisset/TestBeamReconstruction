import numpy as np
import pandas as pd

dataframes = ( pd.read_csv('../Clue/data/output/cluster_data_20_5_20_20_0.csv'),
               pd.read_csv('../Clue/data/output/cluster_data_20_5_20_20_1.csv') )

def plot_energy():
    for df in dataframes:
        mask_clusters = (df['clusterId'] != -1)
        mask_outliers = (df['clusterId'] == -1)
        clusters = df.loc[mask_clusters]
        outliers = df.loc[mask_outliers]

        en_tot = df.loc[:,'weight'].sum()
        en_clusters = clusters.loc[:,'weight'].sum()
        en_outliers = outliers.loc[:,'weight'].sum()
        en_outliers_fraction = en_outliers / en_tot
        print(en_tot, en_clusters, en_outliers, en_outliers_fraction)

        nhits = clusters.groupby('clusterId').size()
        print(nhits)

def plot_position():
    for df in dataframes:
        pass
    

def main():
    #plot_energy()
    plot_position()

if __name__ == '__main__':
    main()
