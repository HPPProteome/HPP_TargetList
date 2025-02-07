#Libraries
from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import gzip
import subprocess
import sys
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


#GENCODE version
version = 47 #Needs to be changed when version updates

#Needed files
gencode_fasta = "gencode.pc_translations.fa.gz"
uniprot_fasta = 'uniprot.fa'
gene_file = "atlasLink.xlsx"

if not os.path.exists(gene_file):
	print(f"{gene_file} file not found, running update_PE.py")
	try:
		subprocess.run(['python3', 'update_PE.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running clean_entries.py: {e}")


gene_file = pd.read_excel(gene_file)
print("Number of genes:", gene_file.shape)

#Changes data from float64s to Strings
gene_file['gencode_symbol'] = gene_file['gencode_symbol'].astype(str)
gene_file['trans_id'] = gene_file['trans_id'].astype(str)
gene_file['CDS'] = gene_file['CDS'].astype(str)

#Downloads GENCODE FASTA file if not availible
if os.path.exists(gencode_fasta):
        print("GENCODE FASTA File Found")
else:
       
        print("Downloading gencode", gencode_fasta)
        url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.pc_translations.fa.gz"
        output_gz_file = "gencode.pc_translations.fa.gz"
        attempt = 0
        max_attempt = 3
        print("Downloading", url)

        while attempt < max_attempt:
            
            try:
                response = requests.get(url, stream=True, timeout=10)
                response.raise_for_status()
                with open(output_gz_file, 'wb') as f:
                        f.write(response.content)
                print("Downloaded", output_gz_file)
                break
            except requests.exceptions.Timeout:
                attempt += 1
                print(f"Request timed out, trying again {attempt}/{max_attempt}")

            except requests.exceptions.RequestException as e:
                print("File not downloaded, error:", e)
                break
        
        if attempt == max_attempt:
            print("File not downloaded, request timed out too many times.")
            sys.exit("Exiting program")


#Downlaods UniProt FASTA file
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



#Creates dictionary to store data from Gencode FASTA file
ensp_dict = {}
symbols_dict = {}
ensg_dict = {}

#Fills 3 dictionaries from FASTA file
#1 that can be indexed by the ENSP number
#1 that can be indexed with the gene symbil
#1 that can be indexed with the ENSG number


with gzip.open(gencode_fasta, "rt") as unzipped_fasta:
    for record in SeqIO.parse(unzipped_fasta, "fasta"):
        header_parts = record.description.split('|')
        ensp_number = header_parts[0]
        gene_symbol = header_parts[-2]
        cds_length = header_parts[-1]
        ensg_number = header_parts[2].split('.')[0]
        enst_number = header_parts[1]
        sequence = str(record.seq)
        ensp_dict[ensp_number] = {"CDS_length": cds_length, "gene_symbol": gene_symbol, "enst":enst_number, "sequence": sequence}
        if gene_symbol not in symbols_dict:
            symbols_dict[gene_symbol] = {"CDS_lengths": [],"ensts": [],"ensps": [], "gene_symbols": [], "sequences":[]}
    
        symbols_dict[gene_symbol]["CDS_lengths"].append(cds_length)
        symbols_dict[gene_symbol]["ensts"].append(enst_number)
        symbols_dict[gene_symbol]["ensps"].append(ensp_number)
        symbols_dict[gene_symbol]['gene_symbols'].append(gene_symbol)
        symbols_dict[gene_symbol]['sequences'].append(sequence)    

        if ensg_number not in ensg_dict:
            ensg_dict[ensg_number]  = {"CDS_lengths": [],"ensts": [],"ensps": [], "gene_symbols": [], "sequences":[]}

        ensg_dict[ensg_number]["CDS_lengths"].append(cds_length)
        ensg_dict[ensg_number]["ensts"].append(enst_number)
        ensg_dict[ensg_number]["ensps"].append(ensp_number)
        ensg_dict[ensg_number]['gene_symbols'].append(gene_symbol)
        ensg_dict[ensg_number]['sequences'].append(sequence) 

#Parses Uniprot FASTA file
uniprot_dict = {}

#Parses FASTA files into dictionaries
for record in SeqIO.parse(uniprot_fasta, "fasta"):
	header_parts = record.description.split('|')
	id = header_parts[1]
	seq = str(record.seq)
	uniprot_dict[id] = seq

#Initializes colunms for the DataFrame <- Not technically neccisary
gene_file['Hydrophobicity'] = ''
gene_file['PI'] = ''
gene_file['Sequence']=''
gene_file['Difference in lengths'] = ''


sameSeq = 0
collected = 0
noUniProtID = 0
noGencode = 0
noMatch = 0
count = 0
incorrectENSP = 0
strippedKey = 0
missed = 0
noSeq = 0

print("Genes to investigate ENSP numbers further")
#Checks to see if the UniProt and Gencode Sequence are the same and assignes it if it is
for index, row in gene_file.iterrows():
	count += 1
#Gets sequence if UniProt and Gencode Match
	if ensp_dict.get(row['ENSP'], {}).get('sequence', "Hello") == uniprot_dict.get(row['uniprot_id'], "World"):
		gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']] #purposely doesn't use .get() to throw error if seq doesn't exist
		gene_file.at[index, 'CDS'] = ensp_dict[row['ENSP']]["CDS_length"]
		gene_file.at[index, 'gencode_symbol'] = ensp_dict[row['ENSP']]["gene_symbol"]
		sameSeq += 1


#Gets ENSP numbers for longest Sequence in Gencode if no UniProtID
	elif row['uniprot_id'] is None or pd.isna(row['uniprot_id']):
                posSeq = ensg_dict.get(row['gene_id'], {}).get('sequences', [])
                max_len = {"len":0, "index":0}
                for i in posSeq:
                        if len(i) > max_len['len']:
                                max_len['len'] = len(i)
                                max_len['index'] = posSeq.index(i)
                i = max_len['index']
                gene_file.at[index, 'Sequence'] = ensg_dict[row['gene_id']]['sequences'][i]
                gene_file.at[index, 'CDS'] = ensg_dict[row['gene_id']]["CDS_lengths"][i]
                gene_file.at[index, 'gencode_symbol'] = ensg_dict[row['gene_id']]['gene_symbols'][i]
                gene_file.at[index, 'ENST'] = ensg_dict[row['gene_id']]['ensts'][i]
                gene_file.at[index, 'ENSP'] = ensg_dict[row['gene_id']]['ensps'][i]
                gene_file.at[index, 'Difference in lengths'] = -int((ensg_dict[row['gene_id']]["CDS_lengths"][i]))
                noUniProtID += 1

#Takes the UniProt ENSP seq from gencode even if Sequences aren't perfect
	elif row['ENSP'] in ensp_dict:
                gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']] #purposely doesn't use .get() to throw error if seq doesn't exist
                gene_file.at[index, 'CDS'] = ensp_dict[row['ENSP']]["CDS_length"]
                gene_file.at[index, 'gencode_symbol'] = ensp_dict[row['ENSP']]["gene_symbol"]
                gene_file.at[index, 'Difference in lengths'] = len(uniprot_dict[row['uniprot_id']]) - int(ensp_dict[row['ENSP']]["CDS_length"])
                incorrectENSP += 1


#Gets the closest length sequence from Gencode to match to UniProte 
	elif uniprot_dict.get(row['uniprot_id']) not in ensg_dict.get(row['gene_id'], {}).get('sequences', []) and uniprot_dict.get(row['uniprot_id']) != None and ensg_dict.get(row['gene_id'], {}).get('sequences', []):
		bestLenMatch = {'len':10**10, 'index':0}
		posSeq = ensg_dict.get(row['gene_id'], {}).get('sequences', False)

		if posSeq:
			for i in posSeq:
				lenDif = abs(len(uniprot_dict[row['uniprot_id']]) - len(i))
				if lenDif < bestLenMatch['len']:
					bestLenMatch['len'], bestLenMatch['index'] = lenDif, posSeq.index(i)
			i = bestLenMatch['index']
			gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']]
			gene_file.at[index, 'CDS'] = ensg_dict[row['gene_id']]["CDS_lengths"][i]
			gene_file.at[index, 'gencode_symbol'] = ensg_dict[row['gene_id']]['gene_symbols'][i]
			gene_file.at[index, 'ENST'] = ensg_dict[row['gene_id']]['ensts'][i]
			gene_file.at[index, 'ENSP'] = ensg_dict[row['gene_id']]['ensps'][i]
			gene_file.at[index, 'Difference in lengths'] = len(uniprot_dict[row['uniprot_id']]) - int((ensg_dict[row['gene_id']]["CDS_lengths"][i]))		
			noMatch += 1

#Tries and get a gencode ENSP number with the same sequence and same ENSG number in UniProt
	elif (row['ENSP'] is None or pd.isna(row['ENSP'])) and not (row['uniprot_id'] is None or pd.isna(row['uniprot_id'])):
		posSeq = ensg_dict.get(row['gene_id'], {}).get('sequences', False)
		if posSeq:
			for i in range(0, len(posSeq)):
				if posSeq[i] == uniprot_dict.get(row['uniprot_id']):
					gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']]
					gene_file.at[index, 'CDS'] = ensg_dict[row['gene_id']]["CDS_lengths"][i]
					gene_file.at[index, 'gencode_symbol'] = ensg_dict[row['gene_id']]['gene_symbols'][i]
					gene_file.at[index, 'ENST'] = ensg_dict[row['gene_id']]['ensts'][i]
					gene_file.at[index, 'ENSP'] = ensg_dict[row['gene_id']]['ensps'][i]
					collected += 1
					break
		else:
			gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']]
			gene_file.at[index, 'Difference in lengths'] = "N/A"
			noGencode += 1


#UniProt only has a stripped ENSP number: ENSP00000494177 instead of ENSP00000494177.1, tries to make the match. Will take first .x number 
	elif "." not in row['ENSP']:
		version = [f".{i}" for i in range(1, 10)]
		for i in version:
			posKey = row['ENSP'] + i
			if posKey in ensp_dict:
				gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']]
				gene_file.at[index, 'CDS'] = ensp_dict[posKey]["CDS_length"]
				gene_file.at[index, 'gencode_symbol'] = ensp_dict[posKey]['gene_symbol']
				if uniprot_dict[row['uniprot_id']] != ensp_dict[posKey]['sequence']:
					gene_file.at[index, 'Difference in lengths'] = len(uniprot_dict[row['uniprot_id']]) - int((ensp_dict[posKey]["CDS_length"]))
				strippedKey += 1

#Just takes the UniProtSeq and ignores ENSP num. For now it is unclear why these ones don't get caught above
	else:
		print(row['gene_id'])
		gene_file.at[index, 'Sequence'] = uniprot_dict[row['uniprot_id']]
		missed += 1		

	if gene_file.at[index, 'Sequence'] == '':
		noSeq += 1

#Calculates Hydrophobicity and pI

	sequence = ''.join([aa for aa in gene_file.at[index, 'Sequence'] if aa in "ACDEFGHIKLMNPQRSTVWY"]) #Removes U (selenocysteine)
	analysis = ProteinAnalysis(sequence)
	gene_file.at[index, 'Hydrophobicity'] = round(analysis.gravy(), 3)
	gene_file.at[index, 'PI'] = round(analysis.isoelectric_point(), 3)


print("\n -----Stats----- \n")
print("Number of UniProt and Gencode Sequence Matches:", sameSeq) 
print("Number of genes with a ENSP seq to match UniProt:", collected)
print("No UniProt ID:", noUniProtID)
print("No Gencode:", noGencode)
print("Incomplete matches:", noMatch)
print("UniProt ENSP not in GENCODE:", incorrectENSP)
print("Number of ENSP numbers with no version:", strippedKey)
print("Genes with just a UniProtSeq (investigate further):", missed)
print("Genes with no sequence:", noSeq)

print("Total genes:", (missed + sameSeq + collected + noUniProtID + noGencode + noMatch + incorrectENSP + strippedKey))

print("\nTrue number of genes:", count)





#Counts and stores outliers to investigate later
noCDS = []
no_uniprot = []
PE5 = []
unreviewed = 0
for index, row in gene_file.iterrows():
        if str(row['isoform']).startswith("ENSG"):
            gene_file.at[index, 'isoform'] = None
        if row['CDS'] == "nan":
            noCDS.append(row)
            gene_file.at[index, 'sequence'] = 'MJA' #Indicates a lack of sequence
        if pd.isnull(row['uniprot_id']):
            no_uniprot.append(row)
        if row['evidence'] == 5:
            PE5.append(row)
        if row['entry_type'] == False:
            unreviewed += 1

print("\nNumber of GENCODE genes not in FASTA", len(noCDS))

print("Number of GENCODE genes with no UniProt connection", len(no_uniprot))
print("Number of unreviewed entries", unreviewed)



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





#Columns that go in the final data frame
columns_to_export = [
    'gene_id', 'gene_name', 'chrom', 'start', 'end', 'transl_type',
    'CDS', 'ENSP', 'ENST',
    'uniprot_id', 'entry_type', 'Reviewed Entry Available', 'entry_name', 'gene_symbol', 'description',
    'protein length', 'evidence', 'Suggested PE', 'Highest nTPM Score','Tissues with nTPM Score Above 1 (/50)', 'found_with', 'isoform', 'refSeq Number', 'Difference in lengths', 
    'EC Number', 'Num Transmembrane Regions', 'Signal Peptide', 'PeptideAtlas Category', 'Observed', 'Distinct', 'Uniquely Mapping', 'Hydrophobicity', 'PI']

#Building dataframes
gene_file_selected = gene_file[columns_to_export]

gene_file_selected = gene_file_selected.rename(columns={"gene_name":"Gene Symbol", "gene_id":"Gene ID", "transl_type": "Translation Type", "gene_symbol":"Uniprot Symbol", "uniprot_id":"UniProtKB ID", "evidence":"PE", "CDS":"CDS Length", "found_with":"Link Made Through", "entry_name":"Entry Name", "chrom":"Chromosome", "description":"Description", "isoform":"Cannonical Isoform","entry_type":"Reviewed"}) 

gene_file_selected['Reviewed'] = gene_file_selected['Reviewed'].astype(bool)

gene_file_selected['Num Transmembrane Regions'] = pd.to_numeric(gene_file_selected['Num Transmembrane Regions'], errors='coerce')
gene_file_selected['Num Transmembrane Regions'] = gene_file_selected['Num Transmembrane Regions'].fillna(0).astype(int).replace(0, pd.NA)
gene_file_selected['Signal Peptide'] = gene_file_selected['Signal Peptide'].replace("['None']", pd.NA) 

#Creates Exception files

#No uniprot entries
noUniprot_df = pd.DataFrame(no_uniprot, columns=columns_to_export)

#Same UniProtID

duplicates = gene_file_selected[gene_file_selected['UniProtKB ID'].notna() & gene_file_selected.duplicated('UniProtKB ID', keep=False)].sort_values(by='UniProtKB ID').reset_index(drop=True)
print("Number of duplicate UniProtKB IDs", duplicates.shape[0])

#Not in Fasta
noCDS_df = pd.DataFrame(noCDS, columns=columns_to_export)


#PE5 Protiens
PE5_df = pd.DataFrame(PE5, columns=columns_to_export)
print("Number of PE 5 proteins:", PE5_df.shape[0])
print("Making Frames")

#Creates table to hold sequences
print("Main table: Supplemental_table_1.xlsx")
output_file = "Supplemental_table_1.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    gene_file_selected.to_excel(writer, sheet_name="Gene Table", index=False)
    secound_sheet_df.to_excel(writer, sheet_name="Column Annotations", index=False)



#Extra tables
noCDS_df.to_excel("No_Fasta_entry.xlsx", index=False)
duplicates.to_excel("Identical_UniProtKB_entries.xlsx",index=False)
PE5_df.to_excel("PE5_Protiens.xlsx", index=False)

print("Tables made")

#Build new FASTA file
identifier_dict = {}
seq_dict = {}

#Builds a fasta file by using the Uniprot ID as an idenfitier.
#A dictionary is used to prevent repeat indentifiers, a ENSP or ENSG number is used in cases of repeats.
#Sequences are assigned to ENSG numbers this way. Uniprot is used when a UniProtKB ID is used, GENCODE is used when a ENSG number is used


def to_fasta(row):


    if row['uniprot_id'] not in identifier_dict and pd.notna(row['uniprot_id']):
        line = f">{row['uniprot_id']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}|{row['description']}|{row['uniprot_id']}|{row['entry_name']}|{row['gene_name']}\n{row['Sequence']}\n"
        identifier_dict[row['uniprot_id']] = "Used"

    elif pd.isna(row['uniprot_id']):
        line = f">{row['ENSP']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}||||{row['gene_name']}\n{row['Sequence']}\n"

    elif pd.isna(row['ENSP']) and row['uniprot_id'] in identifier_dict:
        line = f">{row['gene_id']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}||||{row['gene_name']}\n{row['Sequence']}\n"

    else:
        line = f">{row['ENSP']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}|{row['description']}|{row['uniprot_id']}|{row['entry_name']}|{row['gene_name']}\n{row['Sequence']}\n"

    line = line.replace('nan|','|')
    return line


gene_df = gene_file.sort_values(by='ENSP', na_position='first')

print("Writing FASTA file")
with open('coding_genes.fasta', 'w') as f:
    f.write(''.join(gene_df.apply(to_fasta, axis=1)))
print("File written as coding_genes.fasta")
