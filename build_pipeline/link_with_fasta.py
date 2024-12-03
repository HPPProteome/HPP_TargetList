from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import gzip
import subprocess
import sys

version = 47 #Needs to be changed when version updates

fasta_file = "gencode.pc_translations.fa"
gene_file = "atlasLink.xlsx"

if not os.path.exists(gene_file):
	print(f"{gene_file} file not found, running update_PE.py")
	try:
		subprocess.run(['python3', 'update_PE.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running clean_entries.py: {e}")

gene_file = pd.read_excel(gene_file)

gene_file['gencode_symbol'] = gene_file['gencode_symbol'].astype(str)
gene_file['trans_id'] = gene_file['trans_id'].astype(str)
gene_file['CDS'] = gene_file['CDS'].astype(str)

if os.path.exists(fasta_file):
        print("GENCODE FASTA File Found")
else:
       
        print("Downloading gencode", fasta_file)
        url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.pc_translations.fa.gz"
        output_gz_file = "gencode.pc_translations.gtf.gz"
        output_gtf_file = "gencode.pc_translations.fa"
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

        else:
            print("Unzipping", output_gz_file, "to", output_gtf_file)
            with gzip.open(output_gz_file, 'rb') as f_in:
                     with open(output_gtf_file, 'wb') as f_out:
                             shutil.copyfileobj(f_in, f_out)
            print("Unzipped")


genes_dict = {}
symbols_dict = {}
ensg_dict = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    header_parts = record.description.split('|')
    ensp_number = header_parts[0]
    gene_symbol = header_parts[-2]
    cds_length = header_parts[-1]
    ensg_number = header_parts[2].split('.')[0]
    enst_number = header_parts[1]
    sequence = str(record.seq)
    genes_dict[ensp_number] = {"CDS_length": cds_length, "gene_symbol": gene_symbol, "enst":enst_number, "sequences": sequence}
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
  
gene_file['sequence'] = ''
for index, row in gene_file.iterrows():
	if row['ENSP'] in genes_dict:
		gene_file.at[index, 'CDS'] = genes_dict[row['ENSP']]["CDS_length"]
		gene_file.at[index, 'gencode_symbol'] = genes_dict[row['ENSP']]["gene_symbol"]
		gene_file.at[index, 'sequence'] = genes_dict[row['ENSP']]["sequences"]


for index, row in gene_file.iterrows():
	if row['CDS'] == "nan":
		if row['gene_symbol'] in symbols_dict:
			for i in range(len(symbols_dict[row['gene_symbol']]["ensps"])):
				if row['CDS'] == symbols_dict[row['gene_symbol']]["CDS_lengths"][i]:
					gene_file.at[index, 'CDS'] = symbols_dict[row['gene_symbol']]["CDS_lengths"][i]
					gene_file.at[index, 'gencode_symbol'] = symbols_dict[row['gene_symbol']]['gene_symbols'][i]
					gene_file.at[index, 'ENST'] = symbols_dict[row['gene_symbol']]['ensts'][i]
					gene_file.at[index, 'ENSP'] = symbols_dict[row['gene_symbol']]['ensps'][i]
					gene_file.at[index, 'sequence'] = symbols_dict[row['gene_symbol']]['sequences'][i]

		if row['CDS'] == "nan" and row['gene_symbol'] in symbols_dict:
			list_num = symbols_dict[row['gene_symbol']]["CDS_lengths"]
			best = int(row['protein length'])

			result = [int(x) - best for x in list_num] 
			best_fit = result.index(min(result, key=abs)) 
			
			gene_file.at[index, 'CDS'] = symbols_dict[row['gene_symbol']]["CDS_lengths"][best_fit]
			gene_file.at[index, 'gencode_symbol'] = symbols_dict[row['gene_symbol']]['gene_symbols'][best_fit]
			gene_file.at[index, 'ENST'] = symbols_dict[row['gene_symbol']]['ensts'][best_fit]
			gene_file.at[index, 'ENSP'] = symbols_dict[row['gene_symbol']]['ensps'][best_fit]
			gene_file.at[index, 'sequence'] = symbols_dict[row['gene_symbol']]['sequences'][best_fit]			

for index, row in gene_file.iterrows(): #Takes highest CDS length for any genes that don't have a UniProt Counterpart
	if row['CDS'] == "nan" and  row['gene_id'] in ensg_dict:
		list_num = ensg_dict[row['gene_id']]['CDS_lengths']
		highest_index = list_num.index(max(list_num))
		
		gene_file.at[index, 'CDS'] = ensg_dict[row['gene_id']]['CDS_lengths'][highest_index]
		gene_file.at[index, 'gencode_symbol'] = ensg_dict[row['gene_id']]['gene_symbols'][highest_index]
		gene_file.at[index, 'ENST'] = ensg_dict[row['gene_id']]['ensts'][highest_index]
		gene_file.at[index, 'ENSP'] = ensg_dict[row['gene_id']]['ensps'][highest_index]
		gene_file.at[index, 'sequence'] = ensg_dict[row['gene_id']]['sequences'][highest_index]
	

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

gene_file["CDS"] = pd.to_numeric(gene_file["CDS"], errors="coerce")
gene_file["protein length"] = pd.to_numeric(gene_file["protein length"], errors="coerce")

gene_file["Difference in lengths"] = abs(gene_file["CDS"] - gene_file["protein length"]).replace(0, pd.NA)


columns_to_export = [
    'gene_id', 'gene_name', 'chrom', 'start', 'end', 'transl_type',
    'CDS', 'ENSP', 'ENST',
    'uniprot_id', 'entry_type', 'Reviewed Entry Available', 'entry_name', 'gene_symbol', 'description',
    'protein length', 'evidence', 'Suggested PE', 'Highest nTPM Score','Tissues with nTPM Score Above 1 (/50)', 'found_with', 'isoform', 'Difference in lengths', 
    'EC Number', 'Num Transmembrane Regions', 'Signal Peptide', 'PeptideAtlas Category', 'Observed', 'Distinct', 'Uniquely Mapping', 'sequence']

gene_file_selected = gene_file[columns_to_export]

gene_file_selected = gene_file_selected.rename(columns={"gene_name":"Gene Symbol", "gene_id":"Gene ID", "transl_type": "Translation Type", "gene_symbol":"Uniprot Symbol", "uniprot_id":"UniProtKB ID", "evidence":"PE", "CDS":"CDS Length", "found_with":"Link Made Through", "entry_name":"Entry Name", "chrom":"Chromosome", "description":"Description", "isoform":"Isoform","entry_type":"Reviewed"}) 

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

noUniprot_df.to_excel("No_Uniprot_Entry.xlsx", index=False)
gene_file_selected.to_excel("sequence_table.xlsx", index=False)

gene_file_selected.drop('sequence', axis=1, inplace=True)
gene_file_selected.to_excel(f"Supplemental_table_1_v{version}.xlsx", index=False)

noCDS_df.to_excel("No_Fasta_entry.xlsx", index=False)
duplicates.to_excel("Identical_UniProtKB_entries.xlsx",index=False)
PE5_df.to_excel("PE5_Protiens.xlsx", index=False)

