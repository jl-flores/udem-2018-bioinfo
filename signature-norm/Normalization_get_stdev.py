"""Created:  Tuesday July 10, 2018
   Modified: Monday August 20, 2018
   Jorge Luis Flores
   Computes the standard deviation for the normalizing vector. 
   This can be done either from pre-computed frequency vectors, or by computing 
   the frequency vectors for the motifs in-place."""
import sys
sys.path.append('/u/floresj/Pyth_modules/')
import multiprocessing
from datetime import datetime as dt

from Bio.SeqIO import parse
import motifs_vector_3 as mv
import numpy as np
import pandas as pd

# retrieve the population mean
if __name__ == '__main__':
    MRNA_AVG = pd.read_csv('/u/floresj/mRNA_norm/mrna_avg_vector.csv', index_col=0).as_matrix()[0]

def get_stdev_signs(rna, mcff_args='-t 1', mcff_timeout=300, win_size=79, win_skip=1):
    '''Takes as an argument an RNA as a SeqRecord object.
       Returns a vector of the standard-deviation of signatures for all windows (non-averaged)
       and the number of windows used.
       mcff_args     = mcff parameters, defaulted to -t 1
       mcff_timeout  = time to take before mcff exits
       win_size      = length of windows to use
       win_skip      = how many nucleotides to move to the right each time'''
    
    win_count = 0
    stdev_inside_sum = np.zeros(940)
    
    for i in range(0, len(rna.seq), win_skip):
        # exits loop when the passed window is shorter than win_size
        if len(rna.seq[i:i+win_size]) < win_size:
            break
        
        # calculates the inside sum of (x_i - mu)^2 for the standard deviation
        cur_win = rna.seq[ i:i+win_size ]
        dotbs, shapes = mv.dotbs_and_shapes(cur_win, parameters=mcff_args, timeout_in_seconds=mcff_timeout)
        
        signature_vector  = mv.shape60_ncm40_ncmexp500_expseq340(cur_win, dotbs, shapes)
        stdev_inside_sum += ( (signature_vector - MRNA_AVG) ** 2)
        
        win_count += 1
    
    return ( stdev_inside_sum, win_count )

def get_stdev_sums( rna_freq_df ):
    '''Takes as an argument a df of pre-computed motif frequencies for a single RNA.
       Returns a vector of the "runnning sum" for the standard deviation.'''
    stdev_inside_sum = np.zeros(940)

    # iterate over all the calculated windows
    for ind, rna_row in rna_freq_df.iterrows():
        dotbs       = rna_row.values[0]
        freq_vector = rna_row.values[1:]
        
        # calculates the inside sum of (x_i - mu)^2 for the standard deviation
        signature_vector = freq_vector / dotbs
        stdev_inside_sum += ( (signature_vector - MRNA_AVG) ** 2)
    
    return stdev_inside_sum

if __name__ == '__main__':
    class TimeCounter:
        '''Object used for multithreading'''
        def __init__(self, target):
            self.target = target                # how many TimeCounters must be used
            self.count = 0                      # number of transcripts processed
            self.stdev_sum = np.zeros(940)      # array stores all the results
            #self.nb_windows = 0                 # counts how many windows are used

        def print_progress(self, met_result):
            vector = met_result
            self.stdev_sum += vector            # adds the last result to the sum of all vectors
            #self.nb_windows += nb_win           # keeps track of number of windows used
            self.count += 1                     # keeps track of number of RNAs processed
            
            # status report
            if self.count % 100 == 0 or self.count == self.target:
                print(f'{self.count}\tout of\t{self.target} were processed.')
                
            if self.count % 500 == 0:
                ## keep track of execution time
                with open('/u/floresj/Scripts_log/mrna_normalization_stdev_log.txt', 'a') as fp:
                    fp.write(f'{start_time-dt.now()}\t{self.count}\n')
    
    # initialize logging file
    with open('/u/floresj/Scripts_log/mrna_normalization_stdev_log.txt', 'w') as fp:
        fp.write('Version: Monday August 20, 2018\n')
        fp.write('Time\t\tProcessed\n')
    
    # load transcripts from file
    #all_mrna = []
    #for rna in parse('/u/floresj/Transcriptome_scanning_apta21/GCF_all_mRNA.fa', format='fasta'):
    #    # necessary to transcribe before processing using Mathieu's mv code
    #    rna.seq = rna.seq.transcribe()
    #    all_mrna.append(rna)
    
    ## keep track of execution time
    start_time = dt.now()
    
    # multiprocessing setup
    #multiprocessing.set_start_method('spawn')
    compteur = TimeCounter( 50052 )
    #pool = multiprocessing.Pool(30)
    
    # load all dataframes
    for i in ['01', '02', '03', '04', '05', '06', '07' ,'08', '09', '10', '11', '12', '13']:
        filename = '/u/floresj/mRNA_norm/mRNA_vectors/mrna_folded_int_subset'+i
        freq_df = pd.read_parquet(filename)
        
        # process each df by RNA
        rna_index = list( set(freq_df.index.get_level_values(0)) )
        for rna_name in rna_index:
            #pool.apply_async(get_stdev_sums, ( freq_df.loc[rna_name], ), callback=compteur.print_progress)
            
            compteur.print_progress( get_stdev_sums( freq_df.loc[rna_name] ) )
    
    #pool.close()
    #pool.join()
    
    # load the number of windows used
    nb_windows = 0
    with open('/u/floresj/mRNA_norm/mRNA_vectors/mrna_folded_int_win_dotbs.txt') as fp:
        for line in fp.readlines():
            # skip first three header lines
            if line.startswith('NM'):
                nb_windows += int( line.split('\t')[1].rstrip() )

    # calculate the standard deviation
    vector_stdev = np.sqrt( compteur.stdev_sum / nb_windows )
    print(f'{compteur.count}\t Sequences folded')

    # save the vector in a file, with the shapes as index
    motif_index = mv.list_of_all_motifs()

    vector_df = pd.DataFrame(data=[MRNA_AVG, vector_stdev], columns=motif_index, index=['avg','stdev'])
    vector_df.to_csv('/u/floresj/mRNA_norm/mrna_avg_stdev_vector.csv')