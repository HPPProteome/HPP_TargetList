#Libraries
import pandas as pd
import os
import subprocess
from numpy import nan

#Makes an excel string into a list
def makeList(collection):
	return collection.replace('[', '').replace(']', '').replace("'", '').strip().split(',')


#Table needed 
gene_data = "uniprot_output.xlsx"

#Builds table if table is not availible
if not os.path.exists(gene_data):
	print(f"{gene_data} file not found, running link_to_uniprot.py")
	try:
		subprocess.run(['python3', 'link_to_uniprot.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running link_to_uniprot.py: {e}")

#Reads and modifies gene data for analysis
gene_data = pd.read_excel("uniprot_output.xlsx")
gene_data = gene_data.iloc[:,2:]
gene_data.columns = gene_data.columns.str.strip()

#Since reading an excel file into a dataframe creates a single string, values  need to be re-seperated into lists
columns = ['uniprot_id', 'entry_name', 'gene_symbol', 'description', 'protein length', 'entry_type', 'found_with', 'evidence', 'ENSP', 'ENST', 'isoform', 'Num Transmembrane Regions', 'EC Number', 'Signal Peptide']

for i in columns:
	gene_data[i] = gene_data[i].apply(makeList)

#Optimizes for using gene_id over gene_symbol. Catches any True gene_symbols and False gene_id enteries to store in diff file.
#Will keep False gene_id over True gene_symbol entry.
strange_genes = []

gene_data['Reviewed Entry Available'] = ''
for index, row in gene_data.iterrows():
	found_type = row['found_with']
	contains_id = 'gene_id' in [val.strip() for val in found_type]
	

	if contains_id:
		found_right = [types.strip() == 'gene_id' for types in found_type]
		gene_data.at[index, 'uniprot_id'] = [row['uniprot_id'][i] for i in range(len(found_right)) if found_right[i]]  
		gene_data.at[index, 'entry_name'] = [row['entry_name'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'gene_symbol'] = [row['gene_symbol'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'description'] = [row['description'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'protein length'] = [row['protein length'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'evidence'] = [row['evidence'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'entry_type'] = [row['entry_type'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'ENSP'] = [row['ENSP'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'ENST'] = [row['ENST'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'isoform'] = [row['isoform'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'Num Transmembrane Regions'] = [row['Num Transmembrane Regions'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'EC Number'] = [row['EC Number'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'found_with'] = [row['found_with'][i] for i in range(len(found_right)) if found_right[i]]
		gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i] for i in range(len(found_right)) if found_right[i]]
		for i in range(len(found_type)):
			if row['found_with'][i].strip()  == "gene_name" and row['entry_type'][i].strip():
				strange_genes.append({col: row[col][i] for col in row.index if isinstance(row[col], list)})
				gene_data.at[index, 'Reviewed Entry Available'] = 'yes'

columns = ['uniprot_id', 'entry_name', 'gene_symbol', 'description', 'protein length', 'entry_type', 'evidence', 'found_with', 'isoform', 'EC Number', 'Num Transmembrane Regions']

strange_genes = pd.DataFrame(strange_genes, columns=columns).fillna("")
strange_genes.to_excel("False_ensg_over_name.xlsx", index=False)

#Gets rid of any false entries if a GENCODE gene has a corosponding reviewed UniProt entry
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

        gene_data.at[index, 'Num Transmembrane Regions'] = [row['Num Transmembrane Regions'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        gene_data.at[index, 'EC Number'] = [row['EC Number'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        gene_data.at[index, 'found_with'] = [row['found_with'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
       
#Makes sure that only the lowest  level of exsistance proteins are kept [1,1,4] --> [1,1]
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

		gene_data.at[index, 'EC Number'] = [row['EC Number'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'Num Transmembrane Regions'] = [row['Num Transmembrane Regions'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'found_with'] = [row['found_with'][i].strip() for i in range(len(keeper)) if keeper[i]]
		gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i].strip() for i in range(len(keeper)) if keeper[i]]
#Manually Assigns the correct UniProt ID to 9 protiens
manual_file = "manualFix.tsv"

if not os.path.exists(manual_file):
	print(f"{manual_file} file not found, running manual_fix_file_maker.py")
	try:
		subprocess.run(['python3', 'manual_fix_list_maker.py'], check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error while running manual_fix_list_maker.py: {e}")

manual_fix = pd.read_csv(manual_file, sep='\t')
manual_dict = manual_fix.set_index('Gene ID')['UniProt ID'].to_dict()

for index, row in gene_data.iterrows():
	if row['gene_id'] in manual_dict:
                i = row['uniprot_id'].index(manual_dict[row['gene_id']])
                gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
                gene_data.at[index, 'entry_name'] = row['entry_name'][i]
                gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]
                gene_data.at[index, 'description'] = row['description'][i]
                gene_data.at[index, 'protein length'] = row['protein length'][i]
                gene_data.at[index, 'evidence'] = row['evidence'][i]
                gene_data.at[index, 'entry_type'] = row['entry_type'][i]
                gene_data.at[index, 'ENSP'] = row['ENSP'][i]
                gene_data.at[index, 'ENST'] = row['ENST'][i]
                gene_data.at[index, 'isoform'] = row['isoform'][i]
                gene_data.at[index, 'EC Number'] = row['EC Number'][i]
                gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][i]
                gene_data.at[index, 'found_with'] = row['found_with'][i]
                gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][i]






#Chooses the lowest Isoform if they are present
for index, row in gene_data.iterrows():
	if len(row['isoform']) > 1 and isinstance(row['isoform'], list):
		if not all(item == 'None' or (isinstance(item, str) and item.startswith("ENSG")) for item in row['isoform']):
			useable_isoforms = [[index, isoform] for index, isoform in enumerate(row['isoform']) if isoform != 'None' and not isoform.startswith("ENSG")]
			if not any(None in sublist for sublist in useable_isoforms):
				#print(useable_isoforms)
				#print(row['gene_id'])
				
				i, _ = min(useable_isoforms, key=lambda x: int(x[1].split('-')[-1]))
			
				gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
				gene_data.at[index, 'entry_name'] = row['entry_name'][i]
				gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]  
				gene_data.at[index, 'description'] = row['description'][i]  
				gene_data.at[index, 'protein length'] = row['protein length'][i]
				gene_data.at[index, 'evidence'] = row['evidence'][i]
				gene_data.at[index, 'entry_type'] = row['entry_type'][i]
				gene_data.at[index, 'ENSP'] = row['ENSP'][i]
				gene_data.at[index, 'ENST'] = row['ENST'][i]
				gene_data.at[index, 'isoform'] = row['isoform'][i]
				gene_data.at[index, 'EC Number'] = row['EC Number'][i]
				gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][i]
				gene_data.at[index, 'found_with'] = row['found_with'][i]
				gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][i]



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
                    gene_data.at[index, 'EC Number'] = row['EC Number'][0]
                    gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][0]
                    gene_data.at[index, 'found_with'] = row['found_with'][0]
                    gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][0]
print("Number of repeats:", repeats)

#Chooses highest CDS length when muiltiple entries could be a possible match
for index, row in gene_data.iterrows():
	if len(row['entry_type']) > 1 and not isinstance(row['entry_type'], str):
		mini_dict = {"longest":0, "index":None}
		for i in range(0, len(row['protein length'])):
			if int(row['protein length'][i]) > mini_dict["longest"]:
				mini_dict["longest"] = int(row['protein length'][i])
				mini_dict["index"] = i
	
		i = mini_dict["index"]
		gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
		gene_data.at[index, 'entry_name'] = row['entry_name'][i]
		gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]
		gene_data.at[index, 'description'] = row['description'][i]
		gene_data.at[index, 'protein length'] = row['protein length'][i]
		gene_data.at[index, 'evidence'] = row['evidence'][i]
		gene_data.at[index, 'entry_type'] = row['entry_type'][i]
		gene_data.at[index, 'ENSP'] = row['ENSP'][i]
		gene_data.at[index, 'ENST'] = row['ENST'][i]
		gene_data.at[index, 'isoform'] = row['isoform'][i]
		gene_data.at[index, 'EC Number'] = row['EC Number'][i]
		gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][i]
		gene_data.at[index, 'found_with'] = row['found_with'][i]
		gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][i]

gene_fields = ["gene_id", "gene_name", "chrom", "start", "end", "trans_id", "transl_type", 
    "transl_id", "CDS", "gencode_symbol", "ENSP", "ENST", "uniprot_id", 
    "reviewed", "entry_name", "gene_symbol", "description", "protein length", 
    "entry_type", "evidence", "found_with", "isoform", "EC Number", "Num Transmembrane Regions", "Signal Peptide"]

#Cleans data so that is no longer a list
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
                    gene_data.at[index, 'EC Number'] = row['EC Number'][0]
                    gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][0]
                    gene_data.at[index, 'found_with'] = row['found_with'][0]
                    gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][0]
                    if row['isoform'][0] != None and not row['isoform'][0].startswith("ENSG"):
                        gene_data.at[index, 'isoform'] = row['isoform'][0].strip()
                    else:
                        gene_data.at[index, 'isoform'] = 'None'

print(gene_data.shape)

print("Making Frame: cleaned_table.xlsx")
gene_data.to_excel("cleaned_table.xlsx")
print("done")
