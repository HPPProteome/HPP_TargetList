import pandas as pd


def setOne(collection):
	return collection.replace('[', '').replace(']', '').replace("'", '').strip().split(',')




gene_data = pd.read_excel("uniprot_output.xlsx")
gene_data = gene_data.iloc[:,2:]
gene_data.columns = gene_data.columns.str.strip()

columns = ['uniprot_id', 'entry_name', 'gene_symbol', 'description', 'protein length', 'entry_type', 'evidence', 'ENSP', 'ENST', 'isoform']
for i in columns:
	gene_data[i] = gene_data[i].apply(setOne)

#SETS ENSABL to list

#Gets rid of any false entries in a reviewed set
for index, row in gene_data.iterrows():
    entry_types = row['entry_type']
    contains_true = 'True' in [val.strip() for val in entry_types]

    if contains_true:
        uniprot_ids = row['uniprot_id']
        entry_names = row['entry_name']
        gene_symbols = row['gene_symbol']
        descriptions = row['description']
        protein_lengths = row['protein length']
        evidence = row['evidence']
        entry_types = row['entry_type'] 
        entry_types_bool = [val.strip() == 'True' for val in entry_types]

        filtered_uniprot_ids = [uniprot_ids[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_entry_names = [entry_names[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_gene_symbols = [gene_symbols[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_descriptions = [descriptions[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_protein_lengths = [protein_lengths[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_evidence = [evidence[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        filtered_types = [entry_types[i].strip() for i in range(len(entry_types_bool)) if entry_types_bool[i]]

        gene_data.at[index, 'uniprot_id'] = filtered_uniprot_ids
        gene_data.at[index, 'entry_name'] = filtered_entry_names
        gene_data.at[index, 'gene_symbol'] = filtered_gene_symbols
        gene_data.at[index, 'description'] = filtered_descriptions
        gene_data.at[index, 'protein length'] = filtered_protein_lengths
        gene_data.at[index, 'evidence'] = filtered_evidence
        gene_data.at[index, 'entry_type'] = filtered_types
        
   
        gene_data.at[index, 'ENSP'] = [row['ENSP'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        gene_data.at[index, 'ENST'] = [row['ENST'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        gene_data.at[index, 'isoform'] = [row['isoform'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]

#Makes sure that only the highest level of exsistance protiens are kept
for index, row in gene_data.iterrows():
	if len(row['evidence']) > 1:
		level = [int(exist) for exist in row['evidence']]
		best = min(level)
		keeper = []
		for i in level:
			keeper.append(i==best)
		gene_data.at[index, 'uniprot_id'] = [row['uniprot_id'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'entry_name'] = [row['entry_name'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'gene_symbol'] = [row['gene_symbol'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'description'] = [row['description'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'protein length'] = [row['protein length'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'evidence'] = [row['evidence'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'entry_type'] = [row['entry_type'][i].strip() for i in range(len(keeper)) if keeper[i]]
		
		
		gene_data.at[index, 'ENSP'] = [row['ENSP'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'ENST'] = [row['ENST'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'isoform'] = [row['isoform'][i].strip() for i in range(len(keeper)) if keeper[i]]

repeats = 0
#Uses 1 entry when there are muiltple replicas
for index, row in gene_data.iterrows():
        if isinstance(row['protein length'], list) and len(row['protein length']) > 1:
            stripped_row = [length.strip() for length in row['protein length']]
            if len(set(stripped_row)) == 1:
                    repeats += 1
                    gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][0]
                    gene_data.at[index, 'entry_name'] = row['entry_name'][0]
                    gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][0]
                    gene_data.at[index, 'description'] = row['description'][0]
                    gene_data.at[index, 'protein length'] = row['protein length'][0]
                    gene_data.at[index, 'evidence'] = row['evidence'][0]
                    gene_data.at[index, 'entry_type'] = row['entry_type'][0]

                    gene_data.at[index, 'ENSP'] = row['ENSP'][0]
                    gene_data.at[index, 'ENST'] = row['ENST'][0]
                    gene_data.at[index, 'isoform'] = row['isoform'][0]
print("Number of repeats:", repeats)

#Chooses CDS length when muiltple 
for index, row in gene_data.iterrows():
	if len(row['entry_type']) > 1 and not isinstance(row['entry_type'], str):
		for i in range(0, len(row['protein length'])):
			protein_length = row['protein length'][i].strip()  # Remove extra spaces
			cds_value = str(row['CDS']).split('.')[0].strip()  # Remove extra spaces from CDS and use before '.'
			if protein_length == cds_value:
				#print(len(row['protein length']), len(row['evidence']), type(row['protein length']), type(row['evidence']))
				gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
				gene_data.at[index, 'entry_name'] = row['entry_name'][i]
				gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]
				gene_data.at[index, 'description'] = row['description'][i]
				gene_data.at[index, 'protein length'] = row['protein length'][i]
				gene_data.at[index, 'evidence'] = row['evidence'][i]
				gene_data.at[index, 'entry_type'] = row['entry_type'][i]
				break




gene_fields = ["gene_id", "gene_name", "chrom", "start", "end", "trans_id", "transl_type", 
    "transl_id", "CDS", "gencode_symbol", "ENSP", "ENST", "uniprot_id", 
    "reviewed", "entry_name", "gene_symbol", "description", "protein length", 
    "entry_type", "evidence", "found_with", "isoform"]

#Cleans data that is no longer a list
for index, row in gene_data.iterrows():
	if isinstance(row['entry_name'], list) and len(row['entry_name']) == 1:
                    gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][0]
                    gene_data.at[index, 'entry_name'] = row['entry_name'][0]
                    gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][0]
                    gene_data.at[index, 'description'] = row['description'][0]
                    gene_data.at[index, 'protein length'] = row['protein length'][0]
                    gene_data.at[index, 'evidence'] = row['evidence'][0]
                    gene_data.at[index, 'entry_type'] = row['entry_type'][0]
                    gene_data.at[index, 'ENSP'] = row['ENSP'][0].strip()
                    gene_data.at[index, 'ENST'] = row['ENST'][0].strip()
                    if row['isoform'][0] != None and row['isoform'][0].strip()[:3] != 'ENS':
                        gene_data.at[index, 'isoform'] = row['isoform'][0].strip()
                    else:
                        gene_data.at[index, 'isoform'] = 'None'
hand_look = []
for index, row in gene_data.iterrows():
	if isinstance(row['evidence'],list):
		hand_look.append(row)

human_look = pd.DataFrame(hand_look,columns=gene_fields)
print(human_look.shape)
print(gene_data.shape)

print("Making Frames")
human_look.to_excel("look_over.xlsx")
gene_data.to_excel("full_table.xlsx")
print("done")
