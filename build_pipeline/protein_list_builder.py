from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import sys

version = 47 #Change for formating

gtf_file = "gencode.annotation.gtf.gz"


#Checks for and downloads file
print("Looking for GTF file")
if os.path.exists(gtf_file):
        print("GENCODE GTF File Found")
else:
	print("Downloading", gtf_file, "version", version)
	url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.annotation.gtf.gz"
	output_gz_file = "gencode.annotation.gtf.gz"

	print("Downloading", url)
	max_attempt = 3
	attempt = 0
	while attempt < max_attempt:
		try:
			response = requests.get(url, stream=True, timeout=15)
			response.raise_for_status()

			with open(output_gz_file, 'wb') as f:
				f.write(response.content)
			print("Downloaded", output_gz_file)
			break
		except requests.exceptions.Timeout:
			attempt += 1
			print(f"Attempt to download timed out, retrying {attempt}/{max_attempt}")
		except requests.exceptions.RequestException as e:
			print("An error occured", e) 
			break
	if attempt == max_attempt:
		print(f"Failed to download file after {attempt} attempts")
		sys.exit("Stopping the program.")
	

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
gene_num = 0
gtf_df = gtf_df[gtf_df['feature'].isin(['gene','transcript']) ]
gene_counter = 0
print("Searching for coding genes")
coding_gene = {}
for index, row in gtf_df.iterrows():
	if row["feature"] == 'gene':
		gene_num += 1
		row['attribute'] = parse_attr(row['attribute'])
		gene_type = row['attribute'].get('gene_type', 'N/A')[0]
		if gene_type == 'protein_coding':
			gene_counter += 1
			gene_id = row['attribute'].get('gene_id', 'N/A')[0].split('.')[0]
			gene_name = row['attribute'].get('gene_name', 'N/A')[0]
			chrom = row['chrom'][3:]
			start, end = row['start'], row['end']
			tag = row['attribute'].get('tag', 'N/A')
			trans_type = row['attribute'].get('transcript_type', 'N/A')[0]
			if trans_type not in unexcepted_trans_type:
				if 'readthrough_gene' not in tag:
					coding_gene[gene_id] = {"gene_id": gene_id, "gene_name":gene_name, "chrom": chrom, "start": start, "end": end, "trans_id": None, "transl_type": None, "transl_id": None, "CDS": None}
print(f"There are {gene_num} genes in GENCODE")
print(f"There are {gene_counter} genes with protien coding tag")
print(f"There are {len(coding_gene)} protien coding genes") 
print("Adding translation ids and translation types")
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
final_frame.to_excel("coding_protiens.xlsx") 
print("Done")
