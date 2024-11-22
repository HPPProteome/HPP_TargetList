import os
import requests
import shutil
import gzip
import pandas as pd
import subprocess

#Coverts the protien existance to a numrical figure
def level_converter(description):
	if 'protein level' in description:
		return 1
	elif 'transcript level' in description:
		return 2
	elif 'homology' in description:
		return 3
	elif 'Predicted' in description:
		return 4
	else:
		return 5

def reviewed(checked):
	if checked == 'reviewed':
		return True
	else:
		return False

#Creates a dictionary with the ensamble row, then picks the oldest isoform(if present) and gets the ensg number
def clean_string(s):
    gene_ids = []
    if isinstance(s, str):
        if "\"" in s:
            new_list = s.split("\"")
            for i in range(len(new_list)):
                new_list[i] = new_list[i].split(" ")
            for i in range(len(new_list)):
                for j in range(len(new_list[i])):
                    new_list[i][j] = new_list[i][j].replace(";", "").replace("[", "").replace("]", "")
                    
        else:
            new_list = s.split(' ')
            for j in range(len(new_list)):
                new_list[j] = new_list[j].replace(";", "").replace("[", "").replace("]", "")
        #print(new_list)
        new_list = [sublist for sublist in new_list if any(item != "" for item in sublist)]
        if isinstance(new_list[0], list):
            if len(new_list[0]) > 3:
                for i in range(0, len(new_list)):
                    gene_ids.append({"gene_id":new_list[i][2].split('.')[0], "trans_id":new_list[i][0].split('.')[0], "ensg":new_list[i][1], "isoform": new_list[i][-1]})
        
            else:
                for i in range(0, len(new_list)):
                    gene_ids.append({"gene_id":new_list[i][2].split('.')[0], "trans_id":new_list[i][0].split('.')[0], "ensg":new_list[i][1], "isoform": None})
        elif len(new_list) > 0 and isinstance(new_list[0], str):
            for i in range(0, len(new_list)):
                gene_ids.append({"gene_id":new_list[2].split('.')[0], "trans_id":new_list[0].split('.')[0], "ensg":new_list[1], "isoform": new_list[-1]})
                        
    else:
        gene_ids.append({"gene_id":""})
    return gene_ids

#counts number of transmembrane parts in a protien
def transmem_counter(s):
	if isinstance(s, str):
		return int(s.count("TRANSMEM"))
	else:
		return s



file = "uniprot.tsv"
gene_file = "coding_protiens.xlsx"

if not os.path.exists(gene_file):
	print(f"Missing {gene_file} file, running protein_list_builder.py")
	try:
		subprocess.run(['python3', 'protein_list_builder.py'], check=True)
	except subprocess.CalledProcessError as e:
        	print(f"Error while running protein_list_builder.py: {e}")	

#Checks for and downloads file
print("Looking for uniprot gene file")
if os.path.exists(file):
        print("TSV  File Found")
else:
        print("Downloading gencode.v46.pc_translations.fa")
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Clength%2Cprotein_existence%2Cxref_ensembl_full%2Cid%2Cgene_names%2Cprotein_name%2Cec%2Cft_transmem&format=tsv&query=%28organism_id%3A9606%29"
        output_gz_file = "uni_prot.tsv.gz"
        output_tsv_file = "uniprot.tsv"
        print("Downloading", url)
        response = requests.get(url, stream=True)
        with open(output_gz_file, 'wb') as f:
                f.write(response.content)
        print("Downloaded", output_gz_file)

        print("Unzipping", output_gz_file, "to", output_tsv_file)
        with gzip.open(output_gz_file, 'rb') as f_in:
                with open(output_tsv_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        print("Unzipped")

gene_df = pd.read_excel(gene_file)
id_list = gene_df['gene_id'].tolist()
name_list = gene_df['gene_name'].tolist()
gene_dict = {}
for i in range(len(id_list)):
	gene_dict[id_list[i]] = {"gene_id": id_list[i], "genecode_name": name_list[i], "gencode_symbol": None, "ENSP": [], "ENST":[], "uniprot_id": [], "reviewed": [], "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": [], "isoform":[], "EC Number":[], "Num Transmembrane Regions":[]}
print(len(gene_dict))


all_gene = pd.read_csv(file, sep='\t')
all_gene['Protein existence'] = all_gene['Protein existence'].apply(level_converter)
all_gene['Reviewed'] = all_gene['Reviewed'].apply(reviewed)
all_gene['Transmembrane'] = all_gene['Transmembrane'].apply(transmem_counter)

print("Making connections with gene ids")

count = 0
extra = 0
for index, row in all_gene.iterrows():
	row_dict = clean_string(row['Ensembl'])
	for i in range(0, len(row_dict)):
		gene = row_dict[i]["gene_id"]		
		if gene in gene_dict and row['Reviewed']:
                        gene_dict[gene]['reviewed'].append(row['Reviewed'])
                        gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        gene_dict[gene]['description'].append(row['Protein names']) 
                        gene_dict[gene]['protein length'].append(row['Length'])
                        gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[gene]['found_with'].append("gene_id")
                        gene_dict[gene]['evidence'].append(row['Protein existence'])
                        gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensg"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])

                        gene_dict[gene]['EC Number'].append(row["EC number"])
                        gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])   
                        count +=1
		elif gene in gene_dict and True not in gene_dict[gene]['reviewed']:
                        gene_dict[gene]['reviewed'].append(row['Reviewed'])
                        gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        gene_dict[gene]['description'].append(row['Protein names'])
                        gene_dict[gene]['protein length'].append(row['Length'])
                        gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[gene]['found_with'].append("gene_id")
                        gene_dict[gene]['evidence'].append(row['Protein existence'])
                        gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensg"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])
                        gene_dict[gene]['EC Number'].append(row["EC number"])
                        gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])


print("Making connections with gene symbols")
gene_symbols_dict = {}
for i in gene_dict:
	gene_symbols_dict[gene_dict[i]["genecode_name"]] = i
symbol_count = 0
#Write code to look for genes through symbols; need {"gene_symbol":gene_id}
for index, row in all_gene.iterrows():
	ids = str(row['Gene Names']).replace(';', '').split(' ')
	for name in ids:
		if name.strip() in gene_symbols_dict:
                        ensg = gene_symbols_dict[name.strip()]
                        if gene_dict[ensg]['found_with'] == [] or (True not in gene_dict[ensg]['reviewed'] and row['Reviewed']): 
                                    #print(gene_dict[i]['gencode_symbol'], ids)
                                    gene_dict[ensg]['entry_name'].append(row['Entry Name'])
                                    gene_dict[ensg]['uniprot_id'].append(row['Entry'])
                                    gene_dict[ensg]['description'].append(row['Protein names'])
                                    gene_dict[ensg]['protein length'].append(row['Length'])
                                    gene_dict[ensg]['gene_symbol'].append(name.strip())
                                    gene_dict[ensg]['EC Number'].append(row["EC number"])
                                    gene_dict[ensg]['Num Transmembrane Regions'].append(row["Transmembrane"])
                                    gene_dict[ensg]['found_with'].append("gene_name")
                                    
                                    gene_dict[ensg]['reviewed'].append(row['Reviewed'])
                                    gene_dict[ensg]['evidence'].append(row['Protein existence'])
                                    gene_dict[ensg]['entry_type'].append(row['Reviewed'])

                                    #Code is added for next clean_data.py
                                    gene_dict[ensg]['ENSP'].append(None)
                                    gene_dict[ensg]['ENST'].append(None)
                                    gene_dict[ensg]['isoform'].append(None)
                                    
                                    symbol_count += 1 





not_found = 0

for i in gene_dict:
	if not gene_dict[i]['entry_name']:
		not_found += 1
	if True not in gene_dict[i]["reviewed"]:
		extra += 1

print("Num of reviewed", count)
print("Num of non_reviewed", extra)
#print("Num found with symbol", symbol_count)
print("Num empty", not_found)


print(len(gene_dict))

print("Making Frame")
final_frame = pd.DataFrame(gene_dict).T

gene_fields = ["gene_id", "gene_name", "chrom", "start", "end", "trans_id", "transl_type", 
    "transl_id", "CDS", "gencode_symbol", "ENSP", "ENST", "uniprot_id", 
    "reviewed", "entry_name", "gene_symbol", "description", "protein length", 
    "entry_type", "evidence", "found_with", "isoform"]


merged_df = pd.merge(gene_df, final_frame, on='gene_id', how='inner')

merged_df.to_excel("uniprot_output.xlsx")
