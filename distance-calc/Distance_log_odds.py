"""Created:  Wednesday July 25, 2018
   Modified: Friday July 27, 2018
   Jorge Luis Flores
   Calculates the similarity between each mRNA and a sequence from a given file containing the calculated vectors of mRNA.
   Stores in .csv files a matrix where each row corresponds to an mRNA, and each column corresponds to the sum of log-odds of the 
   79-nucleotide window starting at that nucleotide position (where mRNA is zero-indexed, and each column skips 4 nucleotides).
   
   Similarity is measured by calculating the log odds ratio for each motif, where our test-hypothesis is that the frequency
   of the motif follows a Poisson distribution of lambda=apt21_value, and the null hypothesis is that the frequency follows
   a Poisson distribution of lambda=mrna_avg. This can be written as log( Pr(vector|apt21)/Pr(vector|mrna_avg) )."""
import sys
sys.path.append('/u/floresj/Pyth_modules/')
import multiprocessing
import subprocess

from datetime import datetime as dt

from Bio.SeqIO import parse
import motifs_vector_3 as mv
import numpy as np
import pandas as pd
from math import log
import Distance_measuring as dm

def get_log_odds(vec_i, vec_r, vec_avg, d, n=79, log_base=2):
    """Takes as input two vectors, the vector-of-interest and the reference vector, and calculates the sum of log odds
       of each vector entry as log( Pr( Vi | Vr )/Pr( Vi | Vavg ) ).
       That is, what is the likelihood of the vector-of-interest (Vi) given the reference vector (Vr),
       over the likelihood of the vector-of-interest given the average vector (Vavg).
       
       d = is the number of dotbrackets calculated for the frequencies of vec_i
       n = length of vec_i
       
       Returns a vector of log-odds for each motif."""
    
    dm.same_length(vec_i, vec_r)
    
    # calculate the log-odds for each entry
    log_odds_vec = np.zeros(len(vec_i))
    for i in range(len(vec_i)):
        # get the probability of the vector-of-interest given vec_r, and then given vec_avg
        # the average mRNA vector is in freq/(dotbs*nt), so it must be multiplied by length and dotbs,
        # but the reference vector is only in freq/(dotbs)
        
        if vec_r[i] != 0:
            prob_vec_r   = poisson.pmf(vec_i[i], vec_r[i]*d)
            prob_vec_avg = poisson.pmf(vec_i[i], vec_avg[i]*d*n)
            
            # calculate the log
            #print(f'div\t{(prob_vec_r/prob_vec_avg)}')
            
            # whenever both probabilities are very small, treat as 0/0,
            # so assume that the probability under either model is the same, thus log(odds-ratio)=0
            if (prob_vec_r/prob_vec_avg)==0:
                print(f'Motif\t{ALL_MOTIFS[i]}\n')
                print(f'Frequency\t{vec_i[i]}\t')
                print(f'Pr(vec_i|vec_r)\t{prob_vec_r}')
                print(f'Pr(vec_i|vec_avg)\t{prob_vec_avg}')
                log_odds_vec[i] = 0
            else:
                log_odds_vec[i] = log( (prob_vec_r/prob_vec_avg), log_base )
        # assign as NaN when the motif never occurs in the reference vector
        # done because when lambda=0, Pr(vec_i|vec_r)=0 if the vec_i(val) > 0
        else:
            log_odds_vec[i] = float('NaN')
        
        '''if (prob_vec_r/prob_vec_avg) == 0:
            print(f'Motif\t{ALL_MOTIFS[i]}\n')
            print(f'Frequency\t{vec_i[i]}\t')
            print(f'Pr(vec_i|vec_r)\t{prob_vec_r}')
            print(f'Pr(vec_i|vec_avg)\t{prob_vec_avg}')
            print(f'vec_r\t{vec_r[i]}*{d}\nvec_avg\t{vec_avg[i]}*{d}*{n}')'''
    
    return log_odds_vec

def calculate_similarity(mrna_df, target_signature, rna_id=''):
    """Takes as input the dataframe of a single mRNA, calculates its log-odds ratios, and returns a vector of similarities
       to the passed signature, with the indices*NT_SKIP representing the nucleotide position.
       
       mrna_df          = dataframe of calculated signature (columns) for each nucleotide position (indices). First column is nb_dotbs
       target_signature = p-values of the signature that each window is being compared against."""
    
    
    log_odds_val = []
    # get the log-odds ratio of each row
    for nt_position, freq_vector in mrna_df.iterrows():
        # used to adjust to the number of dotbrackets calculated
        nb_dotbs = freq_vector[0]
        
        cur_logodds = get_log_odds(freq_vector[1:], target_signature, nb_dotbs)
        log_odds_val.append( sum(cur_logodds) )
    
    if rna_id:
        return (log_odds_val, mrna_df.index.values, rna_id)
    else:
        return (log_odds_val, mrna_df.index.values)

if __name__ == '__main__':
    # specify file number to process
    #FILE_NUMBER = input('01, 02, 03, ... 13\nWhich file subset to process...')
    FILE_NUMBER = sys.argv[1]
    print( 'Script is processing subset: '+FILE_NUMBER )
    
    class TimeCounter:
        '''Object used for multithreading'''
        def __init__(self, target):
            self.target = target                # how many TimeCounters must be used
            self.count = 0                      # how many TimeCounters have been used
                                                # in this case, this corresponds to the number of transcripts processed
            self.similarity_dict = {}           # dictionary stores distance per nucleotide, with rna_id as indices

        def print_progress(self, met_result):
            log_odds_val, nt_indices, rna_id = met_result
            
            # add resulting similarity vector to similarity dictionary
            logodds_df = pd.DataFrame(log_odds_val, index=nt_indices).transpose()
            self.similarity_dict[ rna_id ] = logodds_df
            
            self.count += 1                     # keeps track of number of RNAs processed
            
            # prints current status
            if self.count % 100 == 0 or self.count == self.target:
                print(f'{dt.now() - start_time}\t{self.count}')
            # logs into file to keep track
            if self.count % 500 == 0 or self.count == self.target or self.count == 1:
                with open(f'/u/floresj/Scripts_log/log_logodds_file{FILE_NUMBER}.txt', 'w') as fp:
                    fp.write(f'Subset\t{FILE_NUMBER}\n')
                    fp.write(f'Processed {self.count} out of {self.target}\n')
    
    # retrieve the DF of motifs from file
    print('Loading file...')
    df = pd.read_parquet(f'/u/floresj/mRNA_norm/mRNA_vectors/mrna_folded_int_subset{FILE_NUMBER}')
    print('File has loaded.')
    df_indices = sorted( list( {rna_id for rna_id, _ in df.index.values} ) )
    
    # calculate the frequency for aptamer-21
    seq_list = [seq for seq in parse('/u/floresj/Transcriptome_scanning_apta21/aptamer_21.fa', format='fasta')]
    apt_21 = seq_list[0]
    dotbs, shapes = mv.dotbs_and_shapes(apt_21.seq)
    apt21_sign, _, _ = mv.shape60_ncm40_ncmexp500_expseq340(apt_21.seq, dotbs, shapes)
    
    ## keep track of execution time
    print('-------------------------------')
    print('Version: Thursday July 26, 2018')
    print('Time\t\tProcessed')
    start_time = dt.now()
    
    # multiprocessing setup
    multiprocessing.set_start_method('spawn')
    compteur = TimeCounter(len(df_indices))
    pool = multiprocessing.Pool(20)
    
    # process all transcripts
    for ind in df_indices:
        pool.apply_async(calculate_similarity, (df.loc[ind], apt21_sign, ind), callback=compteur.print_progress)
        df.drop(ind, inplace=True)
        #print(ind)
        #compteur.print_progress( calculate_similarity(df.loc[ind], apt21_sign, ind) )
    
    pool.close()
    pool.join()
    
    # save all similarity vectors in a file
    similarities = pd.concat(compteur.similarity_dict)
    
    similarities.to_csv(f'/u/floresj/Transcriptome_scanning_apta21/Distances/logodds_ratio_subset{FILE_NUMBER}.csv')