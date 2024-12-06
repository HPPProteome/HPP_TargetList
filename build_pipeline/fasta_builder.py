import pandas as pd
import os
import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys

#Varable to hold all IDs

identifier_dict = {}

def to_fasta(row):
    if row['UniProtKB ID'] not in identifier_dict and pd.notna(row['UniProtKB ID']):
        line = f">{row['UniProtKB ID']} {row["ENSP"]}|{len(row["sequence"])}|{row['Description']}|{row['UniProtKB ID']}|{row['Entry Name']}|{row['Gene Symbol']}\n{row['sequence']}\n"
        identifier_dict[row['UniProtKB ID']] = "Used"
    elif pd.isna(row['UniProtKB ID']):
        line = f">{row['ENSP']} {row["ENSP"]}|{len(row["sequence"])}||||{row['Gene Symbol']}\n{row["sequence"]}\n"

    elif pd.isna(row['ENSP']) and row['UniProtKB ID'] in identifier_dict:
        line = f">{row['Gene ID']} {row["ENSP"]}|{len(row["sequence"])}||||{row['Gene Symbol']}\n{row["sequence"]}\n"

    else:
        line = f">{row['ENSP']} {row["ENSP"]}|{len(row["sequence"])}|{row['Description']}|{row['UniProtKB ID']}|{row['Entry Name']}|{row['Gene Symbol']}\n{row['sequence']}\n"
    line = line.replace('nan|','|')
    return line

def number(s):
    try:
        int(s)
        return int(s)
    except ValueError:
        if s == "X":
            return 23
        elif s == "Y":
            return 24
        elif s == "M":
            return 25
 
gene_file = "sequence_table.xlsx"
gene_df = pd.read_excel(gene_file)

#Uniprot sequences are taken if Gencode ones do not exist 
uniprot_fasta = 'uniprot.fa'
print("Searching for Uniprot Fasta file")
if not os.path.exists(uniprot_fasta):

	print("Downloading uniprot fasta file for needed sequences")
	url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29"
	attempt = 0
	max_attempt = 3
	
	while attempt < max_attempt:
		try:
			response = requests.get(url, stream=True, timeout=10)
			response.raise_for_status()
			if response.status_code == 200:
				with open('uniprot.fa', 'w') as file:
					file.write(response.text)
				print("FASTA file saved as uniprot.fa")
				break
		except requests.exceptions.Timeout:
			attempt += 1
			print(f"Request to download file timed out, trying again {attempt}/{max_attempt}")
		except requests.exceptions.RequestException as e:
			print("An error occured:", e)
			sys.exit("Exiting program")
else:
	print("File found")
gene_dict = {}

for record in SeqIO.parse(uniprot_fasta, "fasta"):
	header_parts = record.description.split('|')
	id = header_parts[1]
	seq = str(record.seq)
	gene_dict[id] = seq

gene_df['Hydrophobicity'] = ''
gene_df['PI'] = ''

gene_df = gene_df.sort_values(by='UniProtKB ID', na_position='first')

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

gene_df['Chromosome'] = gene_df['Chromosome'].apply(number)
gene_df = gene_df.sort_values(by='Chromosome', ascending=True)
gene_df['Chromosome'] = gene_df['Chromosome'].replace(25, "M").replace(24,"Y").replace(23, "X")

gene_df.drop('sequence', axis=1, inplace=True)
gene_df.to_excel("Supplemental_table_1.xlsx", index=False)
print("Done")
