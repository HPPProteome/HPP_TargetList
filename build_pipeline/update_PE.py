import os
import requests
import zipfile
import pandas as pd

gene_file = "full_table.xlsx"
rna_file = "rna_tissue_consensus.tsv"

print("Looking for RNA_expression file")
if os.path.exists(rna_file):
    print("RNA expressions File Found")
else:
    print("Downloading", rna_file)
    url = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
    output_zip_file = "rna_tissue_consensus.tsv.zip"

    print("Downloading", url)
    response = requests.get(url, stream=True)
    with open(output_zip_file, 'wb') as f:
        f.write(response.content)
    print("Downloaded", output_zip_file)

    print("Unzipping", output_zip_file, "to", rna_file)
    with zipfile.ZipFile(output_zip_file, 'r') as zip_ref:
        zip_ref.extractall(".")
    print("Unzipped")


if not os.path.exists(gene_file):
	print(f"{gene_file} file not found, running clean_entries.p")
	try:
		subprocess.run(['python3', 'clean_entries.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running link_to_uniprot.py: {e}")






rna_data = pd.read_csv(rna_file, sep='\t')
gene_data = pd.read_excel(gene_file)


gene_data['Suggested PE'] = ''

rna_dict = {}

for index, row in rna_data.iterrows():
	gene = row['Gene']
	if gene not in rna_dict:
		rna_dict[gene] = []
	rna_dict[gene].append(row['nTPM'])

count = 0
changed = 0
for index, row in gene_data.iterrows():
	if row['gene_id'] in rna_dict:
		count += 1
		if max(rna_dict[row['gene_id']]) >= 1 and row['evidence'] > 2:
			gene_data.at[index, 'Suggested PE'] = 2
			changed += 1
	#else:
		#print(row['Gene ID'])

print("Making Frame")
gene_data.to_excel("updatedPE.xlsx")
print("done")