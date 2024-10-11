import os
import requests
import shutil
import gzip
import pandas as pd


def level_converter(description):
	if 'protein level' in description:
		return 1
	elif 'transcript level' in description:
		return 2
	elif 'homology' in description:
		return 3
	elif 'predicted' in description:
		return 4
	else:
		return 5

def reviewed(checked):
	if checked == 'reviewed':
		return True
	else:
		return False

def clean_string(s):
    if isinstance(s, str):
        cleaned_str = s.replace(';', ' ').replace('"', ' ').replace('\n', '')
        return cleaned_str
    else:
        return " "
file = "uniprot.tsv"

gene_file = "output.xlsx"

#Checks for and downloads file
print("Looking for uniprot gene file")
if os.path.exists(file):
        print("TSV  File Found")
else:
        print("Downloading gencode.v46.pc_translations.fa")
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Clength%2Cprotein_existence%2Cxref_ensembl_full%2Cid%2Cgene_names%2Cprotein_name&format=tsv&query=%28Human%29+AND+%28model_organism%3A9606%29&sort=protein_name+asc"
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

symbol_list = gene_df['gene_name'].to_list()
print(len(symbol_list))
print(symbol_list[:10])
gene_dict = {}
for i in id_list:
	gene_dict[i] = {"gene_id": i, "gencode_symbol": None, "uniprot_id": [], "reviewed": False, "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": None}
print(len(gene_dict))


for i in range(0, len(gene_dict)):
	gene_dict[id_list[i]]['gencode_symbol'] = symbol_list[i]	


all_gene = pd.read_csv(file, sep='\t')
all_gene['Protein existence'] = all_gene['Protein existence'].apply(level_converter)
all_gene['Reviewed'] = all_gene['Reviewed'].apply(reviewed)
print(all_gene.iloc[80334]['Gene Names'])


print("Making connections with gene ids")

count = 0
extra = 0
for index, row in all_gene.iterrows():
	names = clean_string(row['Ensembl']).split(' ')
	for name in names:
		gene = name.split('.')[0]
		if gene in gene_dict and row['Reviewed']:
                        gene_dict[gene]['reviewed'] = row['Reviewed']
                        gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        gene_dict[gene]['description'].append(row['Protein names']) 
                        gene_dict[gene]['protein length'].append(row['Length'])
                        gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[gene]['found_with'] = "gene_id"
                        gene_dict[gene]['evidence'].append(row['Protein existence'])
                        gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        #print(":)")
                        count += 1
		elif gene in gene_dict and not gene_dict[gene]['reviewed']:
                        gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        gene_dict[gene]['description'].append(row['Protein names'])
                        gene_dict[gene]['protein length'].append(row['Length'])
                        gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[gene]['found_with'] = "gene_id" 
                        gene_dict[gene]['evidence'].append(row['Protein existence'])
                        gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        #print(":(")
print("Making connections with gene symbols")
symbol_count = 0
#Write code to look for genes through symbols
for index, row in all_gene.iterrows():
	ids = str(row['Gene Names'])
	for i in gene_dict:
		if gene_dict[i]['found_with'] == None and gene_dict[i]['gencode_symbol'] in ids:
                        #print(gene_dict[i]['gencode_symbol'], ids)
                        gene_dict[i]['entry_name'].append(row['Entry Name'])
                        gene_dict[i]['uniprot_id'].append(row['Entry'])
                        gene_dict[i]['description'].append(row['Protein names'])
                        gene_dict[i]['protein length'].append(row['Length'])
                        gene_dict[i]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[i]['found_with'] = "gene_name"
                        gene_dict[i]['reviewed'] = row['Reviewed']
                        gene_dict[i]['evidence'].append(row['Protein existence'])
                        gene_dict[i]['entry_type'].append(row['Reviewed'])
                        symbol_count += 1 

not_found = 0
for i in gene_dict:
	if not gene_dict[i]['entry_name']:
		not_found += 1
	if not gene_dict[i]["reviewed"]:
		extra += 1

print("Num of reviewed", count)
print("Num of non_reviewed", extra)
print("Num found with symbol", symbol_count)
print("Num empty", not_found)


print(len(gene_dict))

print("Making Frame")
final_frame = pd.DataFrame(gene_dict).T



merged_df = pd.merge(gene_df, final_frame, on='gene_id', how='inner')

merged_df.to_excel("uniprot_output.xlsx")
