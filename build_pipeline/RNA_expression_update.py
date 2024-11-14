import os
import requests
import zipfile
import pandas as pd

gene_file = "final.xlsx"
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


rna_data = pd.read_csv(rna_file, sep='\t')
gene_data = pd.read_excel(gene_file)

rna_dict = {}

for index, row in rna_data.iterrows():
	gene = row['Gene']
	if gene not in rna_dict:
		rna_dict[gene] = []
	rna_dict[gene].append(row['nTPM'])

count = 0
changed = 0
for index, row in gene_data.iterrows():
	if row['Gene ID'] in rna_dict:
		count += 1
		if max(rna_dict[row['Gene ID']]) >= 1 and row['PE'] > 2:
			gene_data.at[index, 'PE'] = 2 
			changed += 1
	else:
		print(row['Gene ID'])

print("Number of present genes", count)
print("Number of PE scores changed", changed)

print("Making Frame")
gene_data.to_excel("updatedPE.xlsx")
print("done")
