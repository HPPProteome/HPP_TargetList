import HTSeq
import os 
import requests
import pandas as pd
import shutil
import gzip
from Bio import SeqIO

file = "gencode.v46.annotation.gtf" #Change for gencode version"
fasta_file = "gencode.v46.pc_translations.fa"
#Checks for and downloads file
print("Looking for GTF file and FASTA File")
if os.path.exists(file):
        print("GENCODE GTF File Found")
else:
	print("Downloading gencode.v46.pc_translations.fa")
	url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz" #Change URL for new file
	output_gz_file = "gencode.v46.annotation.gtf.gz"
	output_gtf_file = "gencode.v46.annotation.gtf"

	print("Downloading", url)
	response = requests.get(url, stream=True)
	with open(output_gz_file, 'wb') as f:
		f.write(response.content)
	print("Downloaded", output_gz_file)

	print("Unzipping", output_gz_file, "to", output_gtf_file)
	with gzip.open(output_gz_file, 'rb') as f_in:
		with open(output_gtf_file, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	print("Unzipped")

if os.path.exists(fasta_file):
        print("GENCODE FASTA File Found")
else:
        print("Downloading gencode.v46.pc_translations.fa")
        url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.pc_translations.fa.gz" #Change URL for new file
        output_gz_file = "gencode.v46.annotation.gtf.gz"
        output_gtf_file = "gencode.v46.annotation.gtf"

        print("Downloading", url)
        response = requests.get(url, stream=True)
        with open(output_gz_file, 'wb') as f:
                f.write(response.content)
        print("Downloaded", output_gz_file)

        print("Unzipping", output_gz_file, "to", output_gtf_file)
        with gzip.open(output_gz_file, 'rb') as f_in:
                with open(output_gtf_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        print("Unzipped")



valid_chrom = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
print(valid_chrom)
unexcepted_tags = ['readthrough_transcript', 'EnsEMBL_merge_exception', 'readthrough_gene', 'exp_conf', 'RNA_Seq_supported_only']
unexcepted_trans_type = ['nonsense_mediated_decay', 'protein_coding_CDS_not_defined', 'retained_intron', 'TEC']
coding_gene = []
#Gets all coding genes and puts them into a list
cds = 'N/A'
gene_type = None
for feature in read_gtf:
    found = False
    gene_type = feature.attr.get('gene_type', 'N/A')
    gene_id = feature.attr.get('gene_id', 'N/A')
    gene_name = feature.attr.get('gene_name', 'N/A')
    #trans_id = feature.attr.get('transcript_id', 'N/A')
    trans_type = feature.attr.get('transcript_type', 'N/A')
    #transl_id = feature.attr.get('protein_id', 'N/A') 
    chrom = feature.iv.chrom 
    start = feature.iv.start + 1
    end = feature.iv.end 
    tag = feature.attr.get('tag', 'N/A')
    if gene_type == 'protein_coding':
        if chrom in valid_chrom and tag not in unexcepted_tags and trans_type not in unexcepted_trans_type:
            coding_gene.append([gene_id.split('.')[0], gene_name, chrom, start, end, None, trans_type, None, None, tag, gene_type])


print(len(coding_gene))
print(coding_gene[0])
extra = 0
unique_genes = {}
for i in range(0, len(coding_gene)):
	if coding_gene[i][0] in unique_genes:
		extra += 1
	else:
		unique_genes[coding_gene[i][0]] = coding_gene[i]


print("Number of Unique coding genes:", len(unique_genes))
print("Repeats:", extra)

unique_genes_list = []
for i in unique_genes:
	unique_genes_list.append(unique_genes[i])

print("Total genes:", len(unique_genes_list))

fasta_dict = {record.description: record for record in SeqIO.parse(fasta_file, "fasta")}

for i, gene_info in enumerate(unique_genes_list):
    gene_id = gene_info[0]
    for description, data in fasta_dict.items():
        if gene_id in description:
            record_id = data.id.split("|")
            
            unique_genes_list[i][5] = record_id[1].split('.')[0]  # Transcript ID
            unique_genes_list[i][7] = record_id[0].split('.')[0]  # Translation ID
            unique_genes_list[i][8] = record_id[-1]
            break
    
print(unique_genes_list[0])



final_frame = pd.DataFrame(unique_genes_list, columns = ["Gene ID", "Gene Name", "Chromosome", "Start", "End", "Transcript ID", "Transcription Type", "Translation ID", "CDS", "Tag", "Gene Type"])
print(final_frame.head)
final_frame.to_excel("output.xlsx")
