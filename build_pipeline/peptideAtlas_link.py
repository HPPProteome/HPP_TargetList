import requests
import os
import pandas as pd
import subprocess
import sys

url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetCoreProteomeMapping/query_zsun_20241030-091733-695.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_zsun_20241030-091733-695&rs_page_size=1000000&output_mode=tsv"

atlas_file = "peptideAtlas.tsv"
gene_file = "updatedPE.xlsx"

backUp_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins/query_guest_20241206-152057-243.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_guest_20241206-152057-243&rs_page_size=1000000&output_mode=tsv"
backUp_atlas_file = "Secoundary_peptideAtlas.tsv"

attempt = 0
max_attempt = 3
if not os.path.exists(atlas_file):
    while attempt < max_attempt:
        try:
            print("Downloading data from Peptide Atlas")
            response = requests.get(url, stream=True, timeout=10)
            response.raise_for_status()
    
            with open(atlas_file, "wb") as file:
                file.write(response.content)
    
            print(f"File downloaded and saved as {atlas_file}")
            break

        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")

        except requests.exceptions.Timeout:
            attempt += 1
            print(f"Request Timed Out, trying again {attempt}/{max_attempt}")
    if max_attempt == attempt:
        print("Failed to download file")
        sys.exit("Exiting program")

else:
    print("Peptide Atlas file found")
print("Looking for secoundary Atlas File")
attempt = 0
max_attempt = 3
if not os.path.exists(backUp_atlas_file):
    while attempt < max_attempt:
        try:
            print("Downloading data from Peptide Atlas")
            response = requests.get(backUp_url, stream=True, timeout=10)
            response.raise_for_status()

            with open(backUp_atlas_file, "wb") as file:
                file.write(response.content)

            print(f"File downloaded and saved as {backUp_atlas_file}")
            break

        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")

        except requests.exceptions.Timeout:
            attempt += 1
            print(f"Request Timed Out, trying again {attempt}/{max_attempt}")
    if max_attempt == attempt:
        print("Failed to download file")
        sys.exit("Exiting program")

else:
    print("Secoundary Atlas file found")
	


if not os.path.exists(gene_file):
	print(f"Can't find {gene_file}, running update_PE.py")
	try:
		subprocess.run(['python3', 'update_PE.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running update_PE.py: {e}")
print(f"{gene_file} found")

gene_df = pd.read_excel("updatedPE.xlsx")
atlas_df = pd.read_csv(atlas_file, sep='\t')
secound_atlas_df = pd.read_csv(backUp_atlas_file, sep='\t', low_memory=False)

#read into dictionary to go faster
atlas_dict = {}
for index, row in atlas_df.iterrows():
	if isinstance(row['Ensembl_Accession'], str):
		atlas_dict[row['Ensembl_Accession']] = {"catagory":row['PeptideAtlas_Category'], "observed":row['nobs'], "distinct":row['npep'], "unique":row['nunipep'], "uniprot":row['accession']}
secound_atlas_dict = {}
for index, row in secound_atlas_df.iterrows():
	if isinstance(row['biosequence_name'], str):
		secound_atlas_dict[row['biosequence_name']] = {"catagory":row['presence_level'], "observed":row['n_observations'], "distinct":row['n_distinct_peptides']} 



#Set up needed columns
gene_df['PeptideAtlas Category'] = ''
gene_df['Observed'] = ''
gene_df['Distinct'] = ''
gene_df['Uniquely Mapping'] = ''

notIn = 0
mismatch = 0
secound_dict = 0
for index, row in gene_df.iterrows():
	gene = row['gene_id']
	id = row['uniprot_id']
	if gene in atlas_dict:
		gene_df.at[index, 'PeptideAtlas Category'] = atlas_dict[gene]['catagory']
		gene_df.at[index, 'Observed'] = atlas_dict[gene]['observed']
		gene_df.at[index, 'Distinct'] = atlas_dict[gene]['distinct']
		gene_df.at[index, 'Uniquely Mapping'] = atlas_dict[gene]['unique']
		
		if atlas_dict[gene]['uniprot'] != row['uniprot_id']:
			mismatch += 1
	elif id in secound_atlas_dict:
		secound_dict += 1
		gene_df.at[index, 'PeptideAtlas Category'] = secound_atlas_dict[id]['catagory']
		gene_df.at[index, 'Observed'] = secound_atlas_dict[id]['observed']
		gene_df.at[index, 'Distinct'] = secound_atlas_dict[id]['distinct']
	else:
		notIn += 1
	

print(f"There are {mismatch} mismatched uniprot IDs")
print(f"There are {secound_dict} genes found in the secoundary atlas file")

print(f"There are {notIn} genes not in Peptide Atlas")
print("Making Frame")
gene_df.to_excel("atlasLink.xlsx")
print("done")


