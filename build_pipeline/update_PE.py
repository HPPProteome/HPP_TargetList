#Libraries
import os
import requests
import zipfile
import pandas as pd
import subprocess
import sys

#Needed files
gene_file = "cleaned_table.xlsx"
rna_file = "rna_tissue_consensus.tsv.zip"

#Downloads rna_tissue_consensus.tsv from Human Protein Atlas if not present
print("Looking for RNA_expression file")
if os.path.exists(rna_file):
    print("RNA expressions File Found")
else:
    
    print("Downloading", rna_file)
    url = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
    output_zip_file = "rna_tissue_consensus.tsv.zip"

    print("Downloading", url)
    attempt = 0
    max_attempt = 3
    while attempt < max_attempt:
    
        try:
            response = requests.get(url, stream=True, timeout=10)
            with open(output_zip_file, 'wb') as f:
                f.write(response.content)
            print("Downloaded", output_zip_file)
            break
        except requests.exceptions.Timeout:
            attempt += 1
            print(f"Request timed out, trying again {attempt}/{max_attempt}")
        except requests.exceptions.RequestException as e:
            print("Error occured:", e)
            break

    if attempt == max_attempt:
        print("File failed to download")
        sys.exit("Exiting Program")

#Builds cleaned_table.xlsx if not availible 

if not os.path.exists(gene_file):
	print(f"{gene_file} file not found, running clean_entries.py")
	try:
		subprocess.run(['python3', 'clean_entries.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running link_to_uniprot.py: {e}")



#Reads data 
rna_data = pd.read_csv('rna_tissue_consensus.tsv.zip',compression='zip', sep='\t')
gene_data = pd.read_excel(gene_file)

#How high a nTPM score must be to assign a PE2
threshold = 1 #change to be lower or higher

gene_data['Suggested PE'] = ''
gene_data['Highest nTPM Score'] = ''
gene_data[f"Tissues with nTPM Score Above {threshold} (/50)"] = ''
rna_dict = {}

#Collects all nTPM scores for every gene present  
for index, row in rna_data.iterrows():
	gene = row['Gene']
	if gene not in rna_dict:
		rna_dict[gene] = []
	rna_dict[gene].append(row['nTPM'])

count = 0
changed = 0

#Checks to see if a PE > 2 meets the qualifications to be assigned 2.
for index, row in gene_data.iterrows():
	if row['gene_id'] in rna_dict:
		count += 1
		
		gene_data.at[index, 'Highest nTPM Score'] = max(rna_dict[row['gene_id']])
		gene_data.at[index, f"Tissues with nTPM Score Above {threshold} (/50)"]  = len([score for score in rna_dict[row['gene_id']] if score >= threshold])		

		if max(rna_dict[row['gene_id']]) >= threshold and row['evidence'] > 2:
			gene_data.at[index, 'Suggested PE'] = 2
			changed += 1

print("Number of genes in RNA tsv file:", count)
print("Number of genes with new suggested PE score:", changed)

print("Making Frame: updatedPE.xlsx")
gene_data.to_excel("updatedPE.xlsx")
print("done")
