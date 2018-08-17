'''Load all transcripts from the genome file,
   select only for ncRNA that are confirmed
   and store them in a FASTA file.'''
from Bio.SeqIO import parse

genome_filename = '/home/jlflores/Documents/Aptamer21/Genome/GCF_000001405.38_GRCh38.p12_rna.fna'
all_transcripts = [rna for rna in parse(genome_filename, 'fasta')]

# NM_<number> = protein-coding transcript
# NR_<number> = non-protein-coding transcript
# XM_<number> = predicted protein-coding transcript
# XR_<number> = predicted non-protein coding transcript
mrna_count = 0
mrna_predicted_count = 0
ncrna_count = 0
ncrna_predicted_count = 0

access_prefix = {'NM': mrna_count, 'XM': mrna_predicted_count, 'NR': ncrna_count, 'XR': ncrna_predicted_count}
for rna in all_transcripts:
    cur_type = rna.id[:2]
    access_prefix[cur_type] += 1

print(f'{len(all_transcripts)}\tTotal transcripts')
print(f'{access_prefix["NM"]}\tmRNA')
print(f'{access_prefix["XM"]}\tmRNA - Predicted')
print(f'{access_prefix["NR"]}\tncRNA')
print(f'{access_prefix["XR"]}\tncRNA - Predicted')
print('-------------------------------------------------------')

# count types of ncRNAs
ncrna_types = {}
for rna in all_transcripts:
    if rna.id.startswith('NR'):
        # retrieve transcript type and count it
        description = rna.description.split(', ')
        rna_type = description[ len(description)-1 ]
        if rna_type not in ncrna_types:
            ncrna_types[ rna_type ] = 0
        ncrna_types[ rna_type ] += 1

for rna_type in sorted(list(ncrna_types.keys())):
    print(f'{ncrna_types[rna_type]}\t{rna_type}')

# write ncRNAs to file
with open('GCF_all_mRNA.fa', 'w') as fp:
    for rna in all_transcripts:
        if rna.id.startswith('NR'):
            fp.write(f'>{rna.description}\n{str(rna.seq)}\n')