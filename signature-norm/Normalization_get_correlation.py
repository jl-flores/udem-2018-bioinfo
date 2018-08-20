"""Created:  Monday August 13, 2018
   Modified: Monday August 20, 2018
   Jorge Luis Flores
   Generates a normalizing vector of 940 structural motifs from the collection of all mRNAs
   This vector is the average of all windows folded and will be used to normalize"""
import sys
sys.path.append('/u/floresj/Pyth_modules/')

import motifs_vector_3 as mv
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd

from datetime import datetime as dt

ALL_MOTIFS = mv.list_of_all_motifs()

## keep track of execution time
start_time = dt.now()

freq_df = pd.DataFrame()
# load all dataframes
for i in ['01', '02', '03', '04', '05', '06', '07' ,'08', '09', '10', '11', '12', '13']:
    filename = '/u/floresj/mRNA_norm/mRNA_vectors/mrna_folded_int_subset'+i
    tmp_df = pd.read_parquet(filename)
    freq_df = freq_df.append( tmp_df.iloc[:,:101].astype('uint16') )
    
    print(f'{dt.now() - start_time}\tLoaded dataframe {i}')
    
# if all motifs are used, len(freq_df.columns.values) == 941
nb_motifs = len(freq_df.columns.values) - 1

corr_list = []

win_vec = freq_df.iloc[:,0].values

# get all correlation values
for i in range( nb_motifs ):
    # progress report
    print(f'Comparing {ALL_MOTIFS[i]}\tagainst {nb_motifs - (i+1)} other motifs.')
    
    # take the first variable and compare it against other variables, but not to itself
    arr_x = freq_df.iloc[:,i+1].values
    j = i+1
    while j < nb_motifs:
        # retrieve next motif and get correlation
        arr_y = freq_df.iloc[:,j+1].values
        corr_list.append( pearsonr( arr_x/win_vec, arr_y/win_vec ) ) # stores values as tuples
        j += 1
        
    print(f'{dt.now() - start_time}')

# store the condensed correlation/p-val matrix, with the shapes as index
corr_df = pd.DataFrame( data=corr_list, columns=['Pearson r','p-value'] )
corr_df.to_csv(f'/u/floresj/mRNA_norm/mrna_correlation_matrix_condensed_m{nb_motifs}.csv')