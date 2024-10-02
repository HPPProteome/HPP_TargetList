from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import gzip


gtf_file = "gencode.v46.annotation.gtf" #Change for gencode version"
fasta_file = "gencode.v46.pc_translations.fa"


#Checks for and downloads file
print("Looking for GTF file and FASTA File")
if os.path.exists(gtf_file):
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


def parse_attr(attr_str):
	attr = {}
	for attribute in attr_str.split(';'):
		if attribute.strip():  # Ignore empty attributes
			key, value = attribute.strip().split(' ', 1)
			value = value.strip('"')
			if key in attr:
				attr[key].append(value)
			else:
				attr[key] = [value]

	return attr




column_names = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=column_names)
print("GTF file read")

unexcepted_tags = ['readthrough_transcript', 'EnsEMBL_merge_exception', 'readthrough_gene', 'exp_conf', 'RNA_Seq_supported_only']
unexcepted_trans_type = ['nonsense_mediated_decay', 'protein_coding_CDS_not_defined', 'retained_intron', 'TEC']


print("Searching for coding genes")
coding_gene = {}
for index, row in gtf_df.iterrows():
	if row["feature"] == 'gene':
		row['attribute'] = parse_attr(row['attribute'])
		gene_type = row['attribute'].get('gene_type', 'N/A')[0]
		if gene_type == 'protein_coding':
		
			gene_id = row['attribute'].get('gene_id', 'N/A')[0].split('.')[0]
			gene_name = row['attribute'].get('gene_name', 'N/A')[0]
			chrom = row['chrom'][3:]
			start, end = row['start'], row['end']
			tag = row['attribute'].get('tag', 'N/A')
			trans_type = row['attribute'].get('transcript_type', 'N/A')[0]
			if trans_type not in unexcepted_trans_type:
				if 'readthrough_gene' not in tag:
					coding_gene[gene_id] = {"gene_id": gene_id, "gene_name":gene_name, "chrom": chrom, "start": start, "end": end, "trans_id": None, "transl_type": None, "transl_id": None, "CDS": None}
 


fasta_dict = {record.description: record for record in SeqIO.parse(fasta_file, "fasta")}
print(len(coding_gene))
print("Adding trans_id and CDS")
for gene in coding_gene:
	for description, data in fasta_dict.items():
		if gene in description:
			record = data.id.split("|")
			coding_gene[gene]['trans_id'] = record[1].split('.')[0]
			coding_gene[gene]['transl_id'] = record[0].split('.')[0]
			coding_gene[gene]['CDS'] = record[-1]
			break
			

print("Looking for translation type")


#Add transl id, deal with overwritting of tags for diff type so overwrting key/value
for i, row in gtf_df.iterrows():
	translation_type = None
	if row['feature'] == 'transcript':
		row['attribute'] = parse_attr(row['attribute'])
		gene_id = row['attribute'].get('gene_id', 'N/A')[0].split('.')[0]
		if gene_id in coding_gene:
			tag = row['attribute'].get('tag', 'N/A')
			translation_type = None
			
			for i in tag:
				if i == "MANE_Select":
					translation_type = "Mane Select"
				elif i == "Ensembl_canonical":
					translation_type = "canonical"
			if coding_gene[gene_id]['transl_type'] == None:
				coding_gene[gene_id]['transl_type'] = translation_type
	
print("Making Frame")
final_frame = pd.DataFrame(coding_gene).T

print(final_frame.head)
final_frame.to_excel("output.xlsx") 

