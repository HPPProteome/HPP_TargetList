#Libraries
import os
import requests
import shutil
import gzip
import pandas as pd
import subprocess
import re
import sys

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

#Converts reviewed level to a boolean
def reviewed(checked):
	if checked == 'reviewed':
		return True
	else:
		return False

#Creates a dictionary with the ensamble row, then picks the isoform(if present) found in UP000005640_9606.dat and gets the ensp, enst and ensg number
def clean_string(s, id):


	gene_ids = []
	if not isinstance(s, str):
		return gene_ids
	
	elif "[" in s:
		return isoformFinder(s, id)

	else:
		return StringParser(s, id)
	
	
		
def isoformFinder(s, id):
	printed_isoforms = set()
	gene_ids = []
	genes = s.split("\"")
	for i in genes:
		if "[" in i:
			isoform = i.split("[")[1].split("]")[0]
			seperate_ids = i.replace(" ", "").split(";") #ENST, ENSP, ENSG Isoform
			

			if id not in isoformDict:
				if id not in printed_isoforms:
					print(id)
					printed_isoforms.add(id)
				gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":None, "ensp":None, "isoform":seperate_ids[2].split(".")[2].replace("[","").replace("]",""), "refSeq":None})
			elif isoformDict[id] == isoform:
				gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":seperate_ids[0], "ensp":seperate_ids[1], "isoform":isoform, "refSeq":refSeqDict.get(isoform)})	

			else:
				 gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":None, "ensp":None, "isoform":isoformDict.get(id), "refSeq":refSeqDict.get(isoformDict.get(id))})
	return gene_ids

def StringParser(s, id):
	gene_ids = []
	genes = s.split("\"")
	for i in genes:
		if "E" in i: #Checks to make sure the list isn't empty
			seperate_ids = i.replace(" ", "").split(";")
			gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":seperate_ids[0], "ensp":seperate_ids[1], "isoform":isoformDict.get(id), "refSeq":refSeqDict.get(isoformDict.get(id))})
		return gene_ids


#counts number of transmembrane parts in a protien
def transmem_counter(s):
	if isinstance(s, str):
		return int(s.count("TRANSMEM"))
	else:
		return s

#Isolates the n..n part of signal peptides
#SIGNAL 1..22; /evidence="ECO:0000255 --> 1..22
def IsolateSignal(s):
	if not isinstance(s, float):
		num = re.search(r"\d+\.\.\d+", s)
		if num:
			return num.group()
		else:
			return None
	else:
		return None

#Needed files
file = "uniprot.tsv.gz"
gene_file = "coding_protiens.xlsx"
isoformFile = "UP000005640_9606.dat"

#Creates gene_file if it doesn't exist
if not os.path.exists(gene_file):
	print(f"Missing {gene_file} file, running protein_list_builder.py")
	try:
		subprocess.run(['python3', 'protein_list_builder.py'], check=True)
	except subprocess.CalledProcessError as e:
        	print(f"Error while running protein_list_builder.py: {e}")	

#Checks for and downloads UniProtKB  file
print("Looking for uniprot gene file")
if os.path.exists(file):
        print("TSV  File Found")

else:
        print(f"Downloading {file}")
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Clength%2Cprotein_existence%2Cxref_ensembl_full%2Cid%2Cgene_names%2Cprotein_name%2Cft_transmem%2Cec%2Cft_signal&format=tsv&query=%28organism_id%3A9606%29"
        output_gz_file = "uniprot.tsv.gz"
        print("Downloading", url)
        attempt = 0
        max_attempt = 3

        while attempt < max_attempt:
                try:
                        response = requests.get(url, stream=True)
                        response.raise_for_status()
                        with open(output_gz_file, 'wb') as f:
                                f.write(response.content)
                        print("Downloaded", output_gz_file)
                        break
                except requests.exceptions.Timeout:
                        attempt += 1
                        print(f"Request timed out, trying again: {attempt}/{max_attempt}")
                except requests.exceptions.RequestException as e:
                        print("Error", e)
                        break
        if attempt == max_attempt:
                print("Attempt to download file timed out too many times. File not downloaded")
                sys.exit("Exiting Program")



print("Looking for uniprot dat file")

gz_filename = "UP000005640_9606.dat.gz"
output_filename = "UP000005640_9606.dat"

if os.path.exists(output_filename):
    print("DAT File Found")
else:
    print(f"Downloading {gz_filename}")
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.dat.gz"
    print("Downloading", url)
    attempt = 0
    max_attempt = 3

    while attempt < max_attempt:
        try:
            response = requests.get(url, stream=True, timeout=10)
            response.raise_for_status()
            with open(gz_filename, 'wb') as f:
                f.write(response.content)
            print("Downloaded", gz_filename)
            break
        except requests.exceptions.Timeout:
            attempt += 1
            print(f"Request timed out, trying again: {attempt}/{max_attempt}")
        except requests.exceptions.RequestException as e:
            print("Error", e)
            break
    
    if attempt == max_attempt:
        print("Attempt to download file timed out too many times. File not downloaded")
        sys.exit("Exiting Program")
    
    print("Decompressing file")
    with gzip.open(gz_filename, "rb") as f_in:
        with open(output_filename, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed: {output_filename}")






# Creates dataframes and initializes dictionary to hold information
gene_df = pd.read_excel(gene_file)
gene_dict = {}

#Creates a list of all the genes from GENCODE and combines them with their IDs, memory created to hold data collected from UniProtKB
id_list = gene_df['gene_id'].tolist()
name_list = gene_df['gene_name'].tolist()
for i in range(len(id_list)):
	gene_dict[id_list[i]] = {"gene_id": id_list[i], "genecode_name": name_list[i], "gencode_symbol": None, "ENSP": [], "ENST":[], "uniprot_id": [], "reviewed": [], "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": [], "isoform":[], "EC Number":[], "Num Transmembrane Regions":[], "Signal Peptide":[], "refSeq Number":[]}
print(len(gene_dict))

#Reads and modifies UniProtKB data 
all_gene = pd.read_csv(file, sep='\t')

all_gene['Protein existence'] = all_gene['Protein existence'].apply(level_converter)
all_gene['Reviewed'] = all_gene['Reviewed'].apply(reviewed)
all_gene['Transmembrane'] = all_gene['Transmembrane'].apply(transmem_counter)
all_gene['Signal peptide'] = all_gene['Signal peptide'].apply(IsolateSignal)
print("Making connections with gene ids")

#Reads datParser.py file into a dictionary
isoformDict = {} # UniProtID : Cannonical Isoform
refSeqDict = {}  #Isoform : refSequence
with open(isoformFile) as UPfile:
	for line in UPfile:
		if "Sequence=Displayed;" in line:
                        id = line.strip().replace(" ","").replace("CC","").split(";")[0].split("=")[1]
                        uniProtID = id.split("-")[0]
                        isoformDict[uniProtID] = id.split(",")[0]

		if "RefSeq; NP_" in line:
			isoNum = re.search(r"\[(.*?)\]", line)
			if isoNum:
				isoNum = isoNum.group().replace("[","").replace("]","")
				prefix = isoNum.split("-")[0]
				if prefix in isoformDict and isoNum  == isoformDict[prefix]:
					refSeqDict[isoNum] = re.search(r"NP_(.*?);", line).group().replace(";","")
		elif "RefSeq; XP_" in line:
                        isoNum = re.search(r"\[(.*?)\]", line)
                        if isoNum:
                                isoNum = isoNum.group().replace("[","").replace("]","")
                                prefix = isoNum.split("-")[0]
                                if prefix in isoformDict and isoNum  == isoformDict[prefix]:
                                        refSeqDict[isoNum] = re.search(r"XP_(.*?);", line).group().replace(";","")
		elif "RefSeq; NM_" in line:
                        isoNum = re.search(r"\[(.*?)\]", line)
                        if isoNum:
                                isoNum = isoNum.group().replace("[","").replace("]","")
                                prefix = isoNum.split("-")[0]
                                if prefix in isoformDict and isoNum  == isoformDict[prefix]:
                                        refSeqDict[isoNum] = re.search(r"NM_(.*?);", line).group().replace(";","")


#Has correct UniprotKB ID and the ENSG Number, hand chosen gene
exceptions_dict = {"Q6UXT6":"ENSG00000228336"}

#collects data from UniProtKB
count = 0
extra = 0

#Links GENCODE genes to UniProtKB IDs through ENSG numbers 
print("Selected UniProtIDs that have no listed Cannonical Isoform")
for index, row in all_gene.iterrows():
	row_dict = clean_string(row['Ensembl'], row['Entry'])
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
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensp"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])

                        gene_dict[gene]['EC Number'].append(row["EC number"])
                        gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])   
                        gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                        gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"])

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
                        gene_dict[gene]['ENSP'].append(row_dict[i]["ensp"])
                        gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                        gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])
                        gene_dict[gene]['EC Number'].append(row["EC number"])
                        gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])
                        gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                        gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"]) 

	if row['Entry'] in exceptions_dict:
                gene = exceptions_dict[row['Entry']]
                gene_dict[gene]['reviewed'].append(row['Reviewed'])
                gene_dict[gene]['entry_name'].append(row['Entry Name'])
                gene_dict[gene]['uniprot_id'].append(row['Entry'])
                gene_dict[gene]['description'].append(row['Protein names'])
                gene_dict[gene]['protein length'].append(row['Length'])
                gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                gene_dict[gene]['found_with'].append("gene_id")
                gene_dict[gene]['evidence'].append(row['Protein existence'])
                gene_dict[gene]['entry_type'].append(row['Reviewed'])
                gene_dict[gene]['ENSP'].append(None)
                gene_dict[gene]['ENST'].append(None) 
                gene_dict[gene]['isoform'].append(isoformDict.get(row['Entry']))
                gene_dict[gene]['EC Number'].append(row["EC number"])
                gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])
                gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                gene_dict[gene]['refSeq Number'].append(refSeqDict.get(isoformDict.get(row['Entry']))) 


print("Making connections with gene symbols")
gene_symbols_dict = {}
for i in gene_dict:
	gene_symbols_dict[gene_dict[i]["genecode_name"]] = i
symbol_count = 0

#Looks for GENCODE genes' link to UniProtKB through the gene symbols
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
                                    gene_dict[ensg]['isoform'].append(isoformDict.get(row['Entry Name']))
                                    gene_dict[ensg]['Signal Peptide'].append(row["Signal peptide"])
                                    gene_dict[ensg]['refSeq Number'].append(refSeqDict.get(isoformDict.get(row['Entry'])))
                                    symbol_count += 1 




not_found = 0
#Checks for empty ENSG numbers
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

print("Making Frame uniprot_output.xlsx")
final_frame = pd.DataFrame(gene_dict).T

#Gives labels
gene_fields = ["gene_id", "gene_name", "chrom", "start", "end", "trans_id", "transl_type", 
    "transl_id", "CDS", "gencode_symbol", "ENSP", "ENST", "uniprot_id", 
    "reviewed", "entry_name", "gene_symbol", "description", "protein length", 
    "entry_type", "evidence", "found_with", "isoform"]


merged_df = pd.merge(gene_df, final_frame, on='gene_id', how='inner')

merged_df.to_excel("uniprot_output.xlsx")
