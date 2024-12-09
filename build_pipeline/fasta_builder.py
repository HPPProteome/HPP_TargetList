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
Ucount = 0
print("Calculating Hydrophobicity and PI")
#Gets sequence if not present
for index, row in gene_df.iterrows():
	if row['sequence'] == "MJA" and row['UniProtKB ID'] in gene_dict:
		gene_df.at[index, 'sequence'] = gene_dict[row['UniProtKB ID']]

	if row['sequence'] == "MJA" and not row['UniProtKB ID'] in gene_dict:
		print("Not present:", row['Gene ID'], row['Gene Symbol'])
	
	else:
		if 'U' in gene_df.at[index, 'sequence']:
			Ucount += 1

		#Removes U becuase protein analysis can't handle U.
		sequence  = ''.join([aa for aa in gene_df.at[index, 'sequence'] if aa in "ACDEFGHIKLMNPQRSTVWY"])
		#print(sequence)
		analysis = ProteinAnalysis(sequence)
		gene_df.at[index, 'Hydrophobicity'] = round(analysis.gravy(), 3)
		gene_df.at[index, 'PI'] = round(analysis.isoelectric_point(), 3)

print("Number of sequences with U:", Ucount)

print("Writing FASTA file")
with open('coding_genes.fasta', 'w') as f:
    f.write(''.join(gene_df.apply(to_fasta, axis=1)))
print("File written as coding_genes.fasta")


#Includes annotations for what each column means
secound_sheet = {"Gene ID": "ENSG number, from GENCODE GTF.",
    "Gene Symbol": "Gene Symbol, from GENCODE GTF.",
    "Chromosome": "Chromosome location, from GENCODE GTF.",
    "Start": "Nucleotide the gene starts on, from GENCODE FASTA.",
    "End": "Nucleotide the gene ends on, from GENCODE FASTA.",
    "Translation Type": (
        "Represents how a transcript was chosen for a single gene, from GENCODE GTF. "
        "Main Select - Indicates transcript is the single agreed-upon transcript for the protein-coding gene. "
        "Canonical - Indicates it was the 'best' transcript, selected by an individual annotation source."
    ),
    "CDS Length": "Length of gene’s amino acid sequence, from GENCODE FASTA.",
    "ENSP & ENST": "Identifiers taken from UniProtKB. If unavailable, they are taken from GENCODE FASTA.",
    "UniProtKB ID": "Protein Coding Gene UniProtKB ID, taken from UniProtKB.",
    "Reviewed": "Status of entry taken from UniProtKB that matches the gene’s assigned UniProtKB ID.",
    "Reviewed Entry Available": (
        "Indicates that even if a UniProtKB entry with an unreviewed ENSG link was chosen, "
        "a reviewed UniProtKB entry with a gene symbol link is available in False_ensg_over_name.xlsx."
    ),
    "Entry Name": "Gene entry name, from UniProtKB tsv.",
    "Uniprot Symbol": "UniProtKB ID, from UniProtKB tsv.",
    "Description": "Description of UniProtKB entry, from UniProtKB tsv.",
    "Protein Length": "Length of amino acid sequence, from UniProtKB tsv.",
    "PE": (
        "Level of protein existence (1: Evidence at protein level, 2: Evidence at transcript level, "
        "3: Inferred from homology, 4: Protein predicted, and 5: Protein Uncertain), from UniProtKB."
    ),
    "Suggested PE": (
        "Using information from NextProt, if any tissue for a protein coding gene has a score over 1.0, "
        "it suggests a change to 2 if PE score is >2."
    ),
    "Highest nTPM Score": "Provides the highest nTPM score for a gene, from NextProt.",
    "Tissues with nTPM Score Above 1": (
        "Counts the number of tissues with a nTPM score above 1.0. Out of 50 tissues, from NextProt."
    ),
    "Link Made Through": (
        "Indicates whether the UniProtKB entry was linked to the GENCODE entry through the ENG number (gene_id) or its symbol (gene_symbol)."
    ),
    "Isoform": "Notes the specific isoform linked to the ENSG number if present, from UniProtKB tsv.",
    "Difference in Lengths": (
        "Finds the difference in length between the UniProtKB and GENCODE entries’ amino acid sequence."
    ),
    "EC Number": "EC number of the gene, from UniProtKB tsv.",
    "Num Transmembrane Regions": "Number of transmembrane regions, from UniProtKB tsv.",
    "Signal Peptide": "Length of a signal peptide if present, from UniProtKB tsv.",
    "PeptideAtlas Category": "Category of the gene, from PeptideAtlas.",
    "Observed": "Number of observed peptides, from PeptideAtlas.",
    "Distinct": "Number of distinct peptides, from PeptideAtlas.",
    "Uniquely Mapping": "Number of uniquely mapping peptides, from PeptideAtlas.",
    "Hydrophobicity": "Hydrophobicity of amino acid sequence for given gene, calculated with BioPython’s ProteinAnalysis.",
    "PI": "PI number of amino acid sequence for given gene, calculated with BioPython’s ProteinAnalysis."
}




secound_sheet_df = pd.DataFrame(list(secound_sheet.items()), columns=["Column", "Description"])



print("Making updated table")

gene_df['Chromosome'] = gene_df['Chromosome'].apply(number)
gene_df = gene_df.sort_values(by='Chromosome', ascending=True)
gene_df['Chromosome'] = gene_df['Chromosome'].replace(25, "M").replace(24,"Y").replace(23, "X")
gene_df.drop('sequence', axis=1, inplace=True)

output_file = "Supplemental_table_1.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    gene_df.to_excel(writer, sheet_name="Gene Table", index=False)
    secound_sheet_df.to_excel(writer, sheet_name="Column Annotations", index=False)

print("Done")
