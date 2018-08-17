"""Created:  Monday July 30, 2018
   Modified: Monday July 30, 2018
   Jorge Luis Flores
   Calculates the similarity between each mRNA and a sequence from a given file containing the calculated vectors of mRNA.
   Stores in .csv files a matrix where each row corresponds to an mRNA, and each column corresponds to the euclidean distance of the 
   79-nucleotide window starting at that nucleotide position (where mRNA is zero-indexed, and each column skips 4 nucleotides).
   
   Similarity is measured by calculating p-values for each motif, assuming a Poisson distribution. The ln(p-values) are then used
   to calculate an Euclidean distance between the vector of each window and the vector of aptamer-21."""
import sys
sys.path.append('/u/floresj/Pyth_modules/')
import multiprocessing
import subprocess

from datetime import datetime as dt

from Bio.SeqIO import parse
import motifs_vector_3 as mv
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
import Distance_measuring as dm

ALL_MOTIFS = mv.list_of_all_motifs()
WIN_SIZE = 79
#NT_SKIP = 4

if __name__ == '__main__':
    # specify file number to process
    FILE_NUMBER = sys.argv[1]
    print( 'Script is processing subset: '+FILE_NUMBER )

def get_pvalues(mrna_df):
    """Takes the dataframe of a single mRNA, where the rows are the vector containing the frequency of each motif per nucleotide
       (with the index as the nucleotide-position), and the columns are the indices of each motif.
       Returns a list of vectors containing the p-values for each motif per nucleotide position."""
    pval_list = []       # contains the list of signatures as p-values

    # get the p-value of each row while keeping the index of the mRNA
    for nt_position, win_vector in mrna_df.iterrows():
        nb_dotbs = win_vector[0]
        # using survival function instead of 1 - cdf
        cur_poisson, _ = mv.normalize_poisson_precise(win_vector[1:], nb_dotbs, WIN_SIZE)
        pval_list.append( cur_poisson )
    return pval_list

def get_euclidean(mrna_vector, target_sig):
    """Calculates the Euclidean distance of mrna_vector (m) with respect to target_sig (t).
       Each signature is a vector of p-values, and the Euclidean distance is calculated on the ln(p-value)."""
    
    # signatures must be of same length, i.e. 940 entries
    dm.same_length( mrna_vector, target_sig )
    
    # calculate the natural log of each value within the signature array and calculate the Euclidean distance
    log_mrna = np.log(mrna_vector)
    log_targ = np.log(target_sig)
    #for i in range( len(log_mrna) ):
    #    if log_mrna[i] == float('-Inf') or log_mrna[i] == float('Inf'):
    #        print(f'mRNA\t{i}\t{ALL_MOTIFS[i]}\tp-val={mrna_vector[i]}', end='\t|\t')
    #    if log_targ[i] == float('-Inf') or log_targ[i] == float('Inf'):
    #        print(f'Apt21\t{i}\t{ALL_MOTIFS[i]}\tp-val={target_sig[i]}', end='\t|\t')
    
    try:
        euclid_dist = euclidean( np.log(mrna_vector), np.log(target_sig) )
    except ValueError:
        for i in range( len(log_mrna) ):
            if log_mrna[i] == float('-Inf') or log_mrna[i] == float('Inf') or log_targ[i] == float('-Inf') or log_targ[i] == float('Inf'):
                print(f'{log_mrna[i]}')
                print(f'mRNA  {i}\t{ALL_MOTIFS[i]}\tp-val={mrna_vector[i]}')
                print(f'{log_targ[i]}')
                print(f'Apt21 {i}\t{ALL_MOTIFS[i]}\tp-val={target_sig[i]}')
        raise
    
    return euclid_dist
    #return None

def calculate_similarity(mrna_df, target_signature, rna_id=''):
    """Takes as input the dataframe of a single mRNA, calculates its p-values, and returns a vector of distances/similarity
       to the passed signature, with the indices*NT_SKIP representing the nucleotide position.
       
       mrna_df          = dataframe of calculated signature (columns) for each nucleotide position (indices). First column is nb_dotbs
       target_signature = p-values of the signature that each window is being compared against"""
    
    # get p-values for each window within each mRNA
    mrna_pval = get_pvalues(mrna_df)
    
    euclid_values = [ get_euclidean(win_vector, target_signature) for win_vector in mrna_pval ]
    
    if rna_id:
        return (euclid_values, mrna_df.index.values, rna_id)
    else:
        return (euclid_values, mrna_df.index.values)

if __name__ == '__main__':
    class TimeCounter:
        '''Object used for multithreading'''
        def __init__(self, target):
            self.target = target                # how many TimeCounters must be used
            self.count = 0                      # how many TimeCounters have been used
                                                # in this case, this corresponds to the number of transcripts processed
            self.similarity_dict = {}           # dictionary stores distance per nucleotide, with rna_id as indices

        def print_progress(self, met_result):
            euclid_val, nt_indices, rna_id = met_result
            
            # add resulting similarity vector to similarity dictionary
            euclid_val_df = pd.DataFrame(euclid_val, index=nt_indices).transpose()
            self.similarity_dict[ rna_id ] = euclid_val_df
            
            self.count += 1                     # keeps track of number of RNAs processed
            
            # prints current status
            if self.count % 100 == 0 or self.count == self.target:
                print(f'{dt.now() - start_time}\t{self.count}')
            
            # logs into file to keep track
            if self.count % 500 == 0 or self.count == self.target or self.count == 1:
                with open(f'/u/floresj/Scripts_log/log_euclid_dist_{FILE_NUMBER}.txt', 'w') as fp:
                    fp.write(f'Subset\t{FILE_NUMBER}\n')
                    fp.write(f'Processed {self.count} out of {self.target}\n')
    
    # retrieve the DF of motifs from file
    print('Loading file...')
    df = pd.read_parquet(f'/u/floresj/mRNA_norm/mRNA_vectors/mrna_folded_int_subset{FILE_NUMBER}')
    print('File has loaded.')
    df_indices = sorted( list( {rna_id for rna_id, _ in df.index.values} ) )
    
    # calculate the p-values for aptamer-21
    seq_list = [seq for seq in parse('/u/floresj/Transcriptome_scanning_apta21/aptamer_21.fa', format='fasta')]
    apt_21 = seq_list[0]
    dotbs, shapes = mv.dotbs_and_shapes(apt_21.seq)
    a21tmp_signvec, nb_dotbs, _, _ = mv.shape60_ncm40_ncmexp500_expseq340_nodiv(apt_21.seq, dotbs, shapes)
    apt21_signvec, sides = mv.normalize_poisson( a21tmp_signvec, nb_dotbs, len(apt_21.seq) )
    
    ## keep track of execution time
    print('Version: Monday July 30, 2018')
    print('Time\t\tProcessed')
    start_time = dt.now()
    
    # multiprocessing setup
    #multiprocessing.set_start_method('spawn')
    compteur = TimeCounter(len(df_indices))
    #pool = multiprocessing.Pool(20)
    
    # process all transcripts
    for ind in df_indices:
        #pool.apply_async(calculate_similarity, (df.loc[ind], apt21_signvec, ind), callback=compteur.print_progress)
        #df.drop(ind, inplace=True)
        print(ind)
        compteur.print_progress( calculate_similarity(df.loc[ind], apt21_signvec, ind) )
    
    #pool.close()
    #pool.join()
    
    # save all similarity vectors in a file
    distances = pd.concat(compteur.similarity_dict)
    
    distances.to_csv(f'/u/floresj/Transcriptome_scanning_apta21/Distances/euclid_distances_subset{FILE_NUMBER}.csv')