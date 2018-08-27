# udem-2018-bioinfo

Repository for code used during internship at UdeM in summer of 2018. Scripts created by PhD student supervisor are not included. Initial genome dataset (RefSeq version) is available online at https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38.  

# Table of contents
1. [genome-manip](#genome-manip)
2. [signature-norm](#signature-norm)
3. [distance-calc](#distance-calc)

## genome-manip <a name="genome-manip"></a>
Retrieve different subsets of the genome from a local file.
- script1_get_mRNA.py
- script2_get_ncRNA.py

## signature-norm <a name="signature-norm"></a>
Fold fixed-length substructures within transcripts from a local file to calculate to compute baseline statistics for each motif. Save the frequency counts for these substructures into a parquet file.
- Normalization_get_avg.py
- Normalization_fold_all_mrna_ints_upd.py
- Normalization_get_stdev.py
- Normalization_get_correlation.py

## distance-calc <a name="distance-calc"></a>
Calculate the distance between the query RNA and all substructures using different metrics (e.g. Jaccard index, Euclidean distance).
- Distance_measuring.py
- Distance_measuring_pt10_1.py
- Distance_measuring_pt10_2.py
- get_distribution.py
- Distance_euclidean.py
- Distance_log_odds.py
