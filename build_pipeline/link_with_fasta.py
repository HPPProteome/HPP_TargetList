from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import gzip

version = 46

fasta_file = f"gencode.v{version}.pc_translations.fa"
gene_file = pd.read_excel("full_table.xlsx")
gene_file['gencode_symbol'] = gene_file['gencode_symbol'].astype(str)
gene_file['trans_id'] = gene_file['trans_id'].astype(str)
gene_file['CDS'] = gene_file['CDS'].astype(str)

if os.path.exists(fasta_file):
        print("GENCODE FASTA File Found")
else:
        print("Downloading gencode", fasta_file)
        url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.pc_translations.fa.gz"
        output_gz_file = f"gencode.v{version}.pc_translations.gtf.gz"
        output_gtf_file = f"gencode.v{version}.pc_translations.fa"

        print("Downloading", url)
        response = requests.get(url, stream=True)
        with open(output_gz_file, 'wb') as f:
                f.write(response.content)
        print("Downloaded", output_gz_file)

        print("Unzipping", output_gz_file, "to", output_gtf_file)
        with gzip.open(output_gz_file, 'rb') as f_in:
                with open(output_gtf_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        print("Unzipped")


genes_dict = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    header_parts = record.description.split('|')
    ensp_number = header_parts[0]
    gene_symbol = header_parts[-2]
    cds_length = header_parts[-1]
    enst = header_parts[1]
    genes_dict[ensp_number] = {"CDS_length": cds_length, "gene_symbol": gene_symbol, "enst":enst}

symbols_dict = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    header_parts = record.description.split('|')
    ensp = header_parts[0]
    enst = header_parts[1]
    gene_symbol = header_parts[-2]
    cds_length = header_parts[-1]
    if gene_symbol not in symbols_dict:
        symbols_dict[gene_symbol] = {"CDS_lengths": [],"ensts": [],"ensps": [], "gene_symbols": []}
    
    symbols_dict[gene_symbol]["CDS_lengths"].append(cds_length)
    symbols_dict[gene_symbol]["ensts"].append(enst)
    symbols_dict[gene_symbol]["ensps"].append(ensp)
    symbols_dict[gene_symbol]['gene_symbols'].append(gene_symbol)

for index, row in gene_file.iterrows():
	if row['ENSP'] in genes_dict:
		gene_file.at[index, 'CDS'] = genes_dict[row['ENSP']]["CDS_length"]
		gene_file.at[index, 'gencode_symbol'] = genes_dict[row['ENSP']]["gene_symbol"]
for index, row in gene_file.iterrows():
	if row['CDS'] == "nan":
		if row['gene_symbol'] in symbols_dict:
			for i in range(len(symbols_dict[row['gene_symbol']]["ensps"])):
				if row['CDS'] == symbols_dict[row['gene_symbol']]["CDS_lengths"][i]:
					gene_file.at[index, 'CDS'] = symbols_dict[row['gene_symbol']]["CDS_lengths"][i]
					gene_file.at[index, 'gencode_symbol'] = symbols_dict[row['gene_symbol']]['gene_symbols'][i]
					gene_file.at[index, 'ENST'] = symbols_dict[row['gene_symbol']]['ensts'][i]
					gene_file.at[index, 'ENSP'] = symbols_dict[row['gene_symbol']]['ensps'][i]
		
		if row['CDS'] == "nan" and row['gene_symbol'] in symbols_dict:
			list_num = symbols_dict[row['gene_symbol']]["CDS_lengths"]
			best = int(row['protein length'])

			result = [int(x) - best for x in list_num] 
			best_fit = result.index(min(result)) 
			
			gene_file.at[index, 'CDS'] = symbols_dict[row['gene_symbol']]["CDS_lengths"][best_fit]
			gene_file.at[index, 'gencode_symbol'] = symbols_dict[row['gene_symbol']]['gene_symbols'][best_fit]
			gene_file.at[index, 'ENST'] = symbols_dict[row['gene_symbol']]['ensts'][best_fit]
			gene_file.at[index, 'ENSP'] = symbols_dict[row['gene_symbol']]['ensps'][best_fit]
			


num_noCDS = 0
for index, row in gene_file.iterrows():
        if row['CDS'] == "nan":
            num_noCDS += 1
print("Number of no CDS lengths", num_noCDS)

gene_file["CDS"] = pd.to_numeric(gene_file["CDS"], errors="coerce")
gene_file["protein length"] = pd.to_numeric(gene_file["protein length"], errors="coerce")

gene_file["Difference in lengths"] = abs(gene_file["CDS"] - gene_file["protein length"])


columns_to_export = [
    'gene_id', 'gene_name', 'chrom', 'start', 'end', 'transl_type',
    'CDS', 'genecode_name', 'gencode_symbol', 'ENSP', 'ENST',
    'uniprot_id', 'reviewed', 'entry_name', 'gene_symbol', 'description',
    'protein length', 'entry_type', 'evidence', 'found_with', 'isoform', 'Difference in lengths']

gene_file_selected = gene_file[columns_to_export]
gene_file_selected[1:].to_excel("final.xlsx", index=False)




