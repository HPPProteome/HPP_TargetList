import os
import requests
import shutil
import gzip
import pandas as pd

#Coverts the protien existance to a numrical figure
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
                
        oldest = {"num":10*10**10, "index":None}
        if any(isinstance(el, list) for el in new_list):
            for i in range(len(new_list)):
                if len(new_list[i]) > 3:
                    isoform = new_list[i][-1]                
                    if len(isoform) > 0:
                        #print(isoform)
                        if int(isoform[-1]) < oldest["num"]:
                            oldest["num"] = int(new_list[i][-1].split('-')[-1].replace(']', ''))
                            oldest["index"] = i
            if oldest["index"] == None:
                for i in range(0, len(new_list)):
                    if len(new_list[i]) > 1:
                        gene_ids.append({"gene_id":new_list[i][2].split('.')[0], "trans_id":new_list[i][0].split('.')[0], "ensg":new_list[i][1], "isoform": None})
        
            else:
                i = oldest["index"]
                gene_ids.append({"gene_id":new_list[i][2].split('.')[0], "trans_id":new_list[i][0].split('.')[0], "ensg":new_list[i][1], "isoform": new_list[i][-1]})
        else:
            gene_ids.append({"gene_id":new_list[2].split('.')[0], "trans_id":new_list[0].split('.')[0], "ensg":new_list[1], "isoform": new_list[-1]})
    
    else:
        gene_ids.append({"gene_id":""})
    return gene_ids


file = "uniprot.tsv"
gene_file = "coding_protiens.xlsx"

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
name_list = gene_df['gene_name'].tolist()
gene_dict = {}
for i in range(len(id_list)):
	gene_dict[id_list[i]] = {"gene_id": id_list[i], "genecode_name": name_list[i], "gencode_symbol": None, "ENSP": [], "ENST":[], "uniprot_id": [], "reviewed": False, "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": None, "isoform":[]}
print(len(gene_dict))


all_gene = pd.read_csv(file, sep='\t')
all_gene['Protein existence'] = all_gene['Protein existence'].apply(level_converter)
all_gene['Reviewed'] = all_gene['Reviewed'].apply(reviewed)
print(all_gene.iloc[80334])


print("Making connections with gene ids")

count = 0
extra = 0
for index, row in all_gene.iterrows():
	row_dict = clean_string(row['Ensembl'])
	for i in range(0, len(row_dict)):
		gene = row_dict[i]["gene_id"]		
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
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensg"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])
                        count +=1
		elif gene in gene_dict and not gene_dict[gene]['reviewed']:
                        gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        gene_dict[gene]['description'].append(row['Protein names'])
                        gene_dict[gene]['protein length'].append(row['Length'])
                        gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        gene_dict[gene]['found_with'] = "gene_id" 
                        gene_dict[gene]['evidence'].append(row['Protein existence'])
                        gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensg"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])

print("Making connections with gene symbols")
gene_symbols_dict = {}
for i in gene_dict:
	gene_symbols_dict[gene_dict[i]["genecode_name"]] = i
symbol_count = 0
#Write code to look for genes through symbols; need {"gene_symbol":gene_id}
for index, row in all_gene.iterrows():
	ids = str(row['Gene Names']).split(' ')
	for name in ids:
		if name.strip() in gene_symbols_dict:
			ensg = gene_symbols_dict[name.strip()]
			if (not gene_dict[ensg]['reviewed'] and row['Reviewed']) or (gene_dict[ensg]['found_with'] == None): 
                                    #print(gene_dict[i]['gencode_symbol'], ids)
                                    gene_dict[ensg]['entry_name'].append(row['Entry Name'])
                                    gene_dict[ensg]['uniprot_id'].append(row['Entry'])
                                    gene_dict[ensg]['description'].append(row['Protein names'])
                                    gene_dict[ensg]['protein length'].append(row['Length'])
                                    gene_dict[ensg]['gene_symbol'].append(name.strip())
                                    gene_dict[ensg]['found_with'] = "gene_name"
                                    gene_dict[ensg]['reviewed'] = row['Reviewed']
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
	if not gene_dict[i]["reviewed"]:
		extra += 1

print("Num of reviewed", count)
print("Num of non_reviewed", extra)
#print("Num found with symbol", symbol_count)
print("Num empty", not_found)


print(len(gene_dict))

print("Making Frame")
final_frame = pd.DataFrame(gene_dict).T



merged_df = pd.merge(gene_df, final_frame, on='gene_id', how='inner')

merged_df.to_excel("uniprot_output.xlsx")
