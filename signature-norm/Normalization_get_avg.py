"""Created:  Friday July 6, 2018
   Modified: Wednesday July 11, 2018
   Jorge Luis Flores
   Generates a normalizing vector of 940 structural motifs from the collection of all mRNAs
   This vector is the average of all windows folded and will be used to normalize"""
import sys
sys.path.append('/u/floresj/Pyth_modules/')
import multiprocessing
import subprocess

from datetime import datetime as dt
import time
from script1_mp import print_time

from Bio.SeqIO import parse
import motifs_vector_2 as mv
import numpy as np
import pandas as pd


def get_sum_signs(rna, mcff_args='-t 1', mcff_timeout=300, win_size=79, win_skip=1):
    '''Takes as an argument an RNA as a SeqRecord object.
       Returns a vector of the sum of signatures for all windows (non-averaged)
       and the number of windows used.
       mcff_args     = mcff parameters, defaulted to -t 1
       mcff_timeout  = time to take before mcff exits
       win_size      = length of windows to use
       win_skip      = how many nucleotides to move to the right each time'''
    
    # aptamer-21 is 79-nucleotides long
    win_count = 0
    signature_vector = np.zeros(940)
    
    for i in range(0, len(rna.seq), win_skip):
        # exits loop when the passed window is shorter than 79 nt
        if len(rna.seq[i:i+win_size]) < win_size:
            break
        
        # calculates the signature for passed window
        cur_win = rna.seq[i:i+win_size]
        dotbs, shapes = mv.dotbs_and_shapes(cur_win, parameters=mcff_args, timeout_in_seconds=mcff_timeout)
        
        signature_vector += mv.shape60_ncm40_ncmexp500_expseq340(cur_win, dotbs, shapes)
        win_count += 1
    
    return (signature_vector, win_count)

if __name__ == '__main__':
    class TimeCounter:
        '''Object used for multithreading'''
        def __init__(self, target):
            self.starttime = time.time()
            self.target = target                # how many TimeCounters must be used
            self.count = 0                      # how many TimeCounters have been instantiated
                                                # in this case, this corresponds to the number of transcripts processed
            self.vector_sum = np.zeros(940)     # array stores all the results
            self.nb_windows = 0                 # counts how many windows are used

        def print_progress(self, met_result):
            vector, nb_win = met_result
            self.vector_sum += vector           # adds the last result to the sum of all vectors
            self.nb_windows += nb_win           # keeps track of number of windows used
            self.count += 1                     # keeps track of number of RNAs processed
            
            # prints current status
            if self.count % 100 == 0 or self.count == self.target:
                print(f'{self.count}\tout of\t{self.target} were processed.')
                
                ## keep track of execution time
                with open('/u/floresj/mRNA_norm/normalization_fold_all_log.txt', 'a') as fp:
                    fp.write(f'{dt.now() - start_time}\t{self.count}\n')
                
                ## know how much time is left
                #print_time(self.starttime, self.count, self.total)
    
    # load transcripts from file
    all_mrna = []
    for rna in parse('/u/floresj/Transcriptome_scanning_apta21/GCF_all_mRNA.fa', format='fasta'):
        # necessary to transcribe before processing using Mathieu's mv code
        rna.seq = rna.seq.transcribe()
        all_mrna.append(rna)

    # initialize logging file
    with open('/u/floresj/mRNA_norm/normalization_fold_all_log.txt', 'w') as fp:
        fp.write('Version: Wednesday July 11, 2018\n')
        fp.write('Time\t\tProcessed\n')
    
    ## keep track of execution time
    start_time = dt.now()
    
    # multiprocessing setup
    multiprocessing.set_start_method('spawn')
    compteur = TimeCounter(len(all_mrna))
    pool = multiprocessing.Pool(36)

    for rna in all_mrna:
        pool.apply_async(get_sum_signs, (rna,), callback=compteur.print_progress)
        #compteur.print_progress(get_sum_signs(rna))
    pool.close()
    pool.join()
    
    # average the sum of frequencies by number of windows used
    vector_avg = compteur.vector_sum / compteur.nb_windows
    print(f'{compteur.count}\t Sequences folded')

    # save the vector in a file, with the shapes as index
    motif_index = mv.list_of_all_motifs()

    vector_df = pd.DataFrame(vector_avg, index=motif_index).transpose()
    vector_df.to_csv('/u/floresj/mRNA_norm/mrna_avg_vector.csv')
