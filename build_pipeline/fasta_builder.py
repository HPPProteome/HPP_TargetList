import pandas as pd
import os
import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def to_fasta(row):
    return f">{row['UniProtKB ID']} {row["ENSP"]}|{len(row["sequence"])}|{row['Description']}|{row['Entry Name']}|{row['Gene Symbol']}\n{row['sequence']}\n"

gene_file = "sequence_table.xlsx"
gene_df = pd.read_excel(gene_file)

#Uniprot sequences are taken if Gencode ones do not exist 
uniprot_fasta = 'uniprot.fa'
if not os.path.exists(uniprot_fasta):
	print("Downloading uniprot fasta file for needed sequences")
	url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29"
	response = requests.get(url)

	if response.status_code == 200:
		with open('uniprot.fa', 'w') as file:
			file.write(response.text)
		print("FASTA file saved as uniprot.fa")
	else:
		print("Failed to download the file.")

gene_dict = {}

for record in SeqIO.parse(uniprot_fasta, "fasta"):
	header_parts = record.description.split('|')
	id = header_parts[1]
	seq = str(record.seq)
	gene_dict[id] = seq

gene_df['Hydrophobicity'] = ''
gene_df['PI'] = ''



print("Calculating Hydrophobicity and PI")
#Gets sequence if not present
for index, row in gene_df.iterrows():
	if row['sequence'] == "MJA" and row['UniProtKB ID'] in gene_dict:
		gene_df.at[index, 'sequence'] = gene_dict[row['UniProtKB ID']]

	if row['sequence'] == "MJA" and not row['UniProtKB ID'] in gene_dict:
		print("Not present:", row['Gene ID'], row['Gene Symbol'])
	
	else:

		#Removes U becuase protein analysis can't handle U.
		sequence  = ''.join([aa for aa in gene_df.at[index, 'sequence'] if aa in "ACDEFGHIKLMNPQRSTVWY"])
		#print(sequence)
		analysis = ProteinAnalysis(sequence)
		gene_df.at[index, 'Hydrophobicity'] = round(analysis.gravy(), 3)
		gene_df.at[index, 'PI'] = round(analysis.isoelectric_point(), 3)



print("Writing FASTA file")
with open('coding_genes.fasta', 'w') as f:
    f.write(''.join(gene_df.apply(to_fasta, axis=1)))
print("File written as coding_genes.fasta")

print("Making updated table")

gene_df.drop('sequence', axis=1, inplace=True)
gene_df.to_excel("Updated_Sup_table.xlsx")
print("Done")
