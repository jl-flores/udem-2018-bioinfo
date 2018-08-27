"""Created:  Tuesday July 10, 2018
   Modified: Monday August 20, 2018
   Jorge Luis Flores
   Generates 13 dataframes containing vectors for all calculated windows for all mRNAs
   These vectors contain nb_dotbs, [ALL_MOTIFS]"""
import sys
sys.path.append('/u/floresj/Pyth_modules/')
import multiprocessing

from datetime import datetime as dt

from Bio.SeqIO import parse
import motifs_vector_3 as mv
import numpy as np
import pandas as pd

ALL_MOTIFS = mv.list_of_all_motifs()
TRANSCRIPTS_PER_FILE = 4000
WIN_SKIP = 4
MCFF_T = '-t 1'

def get_sum_signs_int(rna, mcff_args='-t 1', mcff_timeout=300, win_size=79, win_skip=WIN_SKIP):
    '''Takes as an argument an RNA as a SeqRecord object.
       Returns a list of vectors for the motif counts, along with the windows used and the number of dot-brackets used.
       Also returns the rna_id of the processed transcript.
       Used to make more lightweight .csv files in which all numeric entries are ints
       
       mcff_args     = mcff parameters, defaulted to -t 1
       mcff_timeout  = time to take before mcff exits
       win_size      = length of windows to use
       win_skip      = how many nucleotides to move to the right each time'''
    
    win_count = 0
    all_sign_vectors = []
    dotbs_in_win = []
    nt_position = []
    
    for i in range(0, len(rna.seq), win_skip):
        # exits loop when the passed window is shorter than win_size
        if len(rna.seq[i:i+win_size]) < win_size:
            break
        
        # calculate the signature for passed window
        cur_win = rna.seq[i:i+win_size]
        dotbs, shapes = mv.dotbs_and_shapes(cur_win, parameters=mcff_args, timeout_in_seconds=mcff_timeout)
        
        cur_signature_vector, nb_dotbs, _, _ = mv.shape60_ncm40_ncmexp500_expseq340_nodiv(cur_win, dotbs, shapes)
        win_count += 1
        
        all_sign_vectors.append(cur_signature_vector)
        
        # keeps track of the dotbs and nt position
        dotbs_in_win.append(nb_dotbs)
        nt_position.append(i)
    
    int_vectors = np.array( [arr.astype(int) for arr in all_sign_vectors] )
    
    return (int_vectors, win_count, dotbs_in_win, rna.id, nt_position)

if __name__ == '__main__':
    class TimeCounter:
        '''Object used for multithreading'''
        def __init__(self, target):
            self.target = target                # how many TimeCounters must be used
            self.count = 0                      # number of transcripts processed
            self.all_df = {}                    # stores the dataframe of each mRNA

        def print_progress(self, met_result):
            vector_list, nb_win, dotbs_in_win, rna_name, nt_position = met_result
            
            # store the vectors in a file
            # add columns for nb_dotbs
            vector_list = np.insert(vector_list, 0, dotbs_in_win, axis=1)
            
            # write nb_win to a .txt file
            with open('/u/floresj/mRNA_norm/mRNA_vectors/updated/mrna_folded_t5_nbwin.txt', 'a') as fp:
                fp.write(f'{rna_name}\t{nb_win}\n')
            
            # use nt_position as index
            vector_df = pd.DataFrame( vector_list, columns = ['dotbs_in_win']+ALL_MOTIFS, index=nt_position, dtype=np.uint16 )
            self.all_df[rna_name] = vector_df
            
            self.count += 1                     # keeps track of number of RNAs processed
            
            # prints current status
            if self.count % 100 == 0 or self.count == self.target:
                print(f'{dt.now() - start_time}\t{self.count}')
            
            # outputs to file once enough mRNAs have been processed
            if ( (self.count % TRANSCRIPTS_PER_FILE == 0) or (self.count == self.target) ):
                total_df = pd.concat( compteur.all_df )
                if self.count % TRANSCRIPTS_PER_FILE == 0:
                    total_df.to_parquet(f'/u/floresj/mRNA_norm/mRNA_vectors/updated/mrna_folded_t1_subset{int(self.count/TRANSCRIPTS_PER_FILE)}')
                else:
                    total_df.to_parquet(f'/u/floresj/mRNA_norm/mRNA_vectors/updated/mrna_folded_t1_subset{int(self.count//TRANSCRIPTS_PER_FILE)+1}')
                
                # flushes dictionary
                total_df = pd.DataFrame()
                self.all_df = {}
                print('\tTranscripts were saved to file, and memory was flushed.')
    
    print('Version: Monday August 20, 2018')
    print('Time\t\tProcessed')
    
    # file to store nb_win
    with open('/u/floresj/mRNA_norm/mRNA_vectors/updated/mrna_folded_t5_nbwin.txt', 'w') as fp:
        fp.write('Version: Monday August 20, 2018\n')
        fp.write(f'Windows calculated every {WIN_SKIP} nucleotides, using {MCFF_T}.\n')
        fp.write('rna_id\tnb_win\n')
    
    # load transcripts from file
    all_mrna = []
    for rna in parse('/u/floresj/Transcriptome_scanning_apta21/step1_Normalization_and_folding/GCF_all_mRNA.fa', format='fasta'):
        # necessary to transcribe before processing using Mathieu's mv code
        rna.seq = rna.seq.transcribe()
        all_mrna.append(rna)
    
    ## keep track of execution time
    start_time = dt.now()
    
    # multiprocessing setup
    compteur = TimeCounter(len(all_mrna))
    multiprocessing.set_start_method('spawn')
    pool = multiprocessing.Pool(40)

    for rna in all_mrna:
        pool.apply_async(get_sum_signs_int, (rna, MCFF_T, 500), callback=compteur.print_progress)
        #compteur.print_progress(get_sum_signs(rna))
    pool.close()
    pool.join()
