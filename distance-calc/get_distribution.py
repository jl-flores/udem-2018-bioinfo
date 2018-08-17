import numpy as np
import pandas as pd

from Bio.SeqIO import parse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# load all dataframes and concatenate into one
df_total = pd.DataFrame()
for i in ['01', '02', '03', '04', '05', '06', '07' ,'08', '09', '10.1', '10.2', '11', '12', '13']:
    filename = '/u/floresj/Transcriptome_scanning_apta21/Distances/similarities_subset'+i+'.csv'
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=1)
    df.drop(1, axis=1, inplace=True)
    
    df_total = df_total.append(df)
    #df_total = df_total.append(df, sort=False)
    
    print(i)
    
from datetime import datetime as dt
start_time = dt.now()
# retrieve all similarity values
val_list = []
count = 0
for rna_id, row in df_total.iterrows():
    # store all the calculated Jaccard index values
    val_list = val_list + list(row.dropna())
    
    count+=1
    if count%5000==0:
        print(str(dt.now() - start_time)+'\t'+str(count))
        
# plot the values in a histogram
for BINS in [10, 50, 100]:
    plt.hist(val_list, bins=BINS)
    
    plt.title('Distribution of Jaccard indices for all windows')
    plt.xlabel('Jaccard index')
    plt.ylabel('Frequency')
    
    plt.savefig(f'/u/floresj/Transcriptome_scanning_apta21/Distances_analysis/jacc_index_distr_pval0.5_bin{BINS}.pdf')
    
    plt.close()