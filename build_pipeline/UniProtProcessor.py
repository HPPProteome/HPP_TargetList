import os
import requests
import shutil
import gzip
import pandas as pd
import subprocess
import re
import sys

class UniProtProcessor:
    def __init__(self):
        self.uniprot_file = "uniprot.tsv.gz"
        self.gene_file = "protien_coding_genes.xlsx"
        self.isoformFile = "Uniprot.dat"
        self.output_file = "uniprot_output.xlsx"
        
        self.gene_df = pd.read_excel(self.gene_file)
        self.uniprot_genes = pd.DataFrame()
        
        self.gene_dict = {}
        self.key_words = {} #UniProtID (s) : Associated tags
        self.isoformDict = {} # UniProtID : Cannonical Isoform
        self.refSeqDict = {}  #Isoform : refSequence


    def level_converter(self, description):
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
        
    def reviewed(self, checked):
        if checked == 'reviewed':
            return True
        else:
            return False
        
    #Creates a dictionary with the ensamble row, then picks the isoform(if present) found in UP000005640_9606.dat and gets the ensp, enst and ensg number
    def clean_string(self, s, id):
        gene_ids = []
        if not isinstance(s, str):
            return gene_ids
        
        elif "[" in s:
            return self.isoformFinder(s, id)

        else:
            return self.StringParser(s, id)
	
	
		
    def isoformFinder(self, s, id):
        printed_isoforms = set()
        gene_ids = []
        genes = s.split("\"")
        for i in genes:
            if "[" in i:
                isoform = i.split("[")[1].split("]")[0]
                seperate_ids = i.replace(" ", "").split(";") #ENST, ENSP, ENSG Isoform
                

                if id not in self.isoformDict:
                    if id not in printed_isoforms:
                        print(id)
                        printed_isoforms.add(id)
                    gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":None, "ensp":None, "isoform":seperate_ids[2].split(".")[2].replace("[","").replace("]",""), "refSeq":None})
                elif self.isoformDict[id] == isoform:
                    gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":seperate_ids[0], "ensp":seperate_ids[1], "isoform":isoform, "refSeq":self.refSeqDict.get(isoform)})	

                else:
                    gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":None, "ensp":None, "isoform":self.isoformDict.get(id), "refSeq":self.refSeqDict.get(self.isoformDict.get(id))})
        return gene_ids

    def StringParser(self, s, id):
        gene_ids = []
        genes = s.split("\"")
        for i in genes:
            if "E" in i: #Checks to make sure the list isn't empty
                seperate_ids = i.replace(" ", "").split(";")
                gene_ids.append({"gene_id":seperate_ids[2].split(".")[0], "trans_id":seperate_ids[0], "ensp":seperate_ids[1], "isoform":self.isoformDict.get(id), "refSeq":self.refSeqDict.get(self.isoformDict.get(id))})
        return gene_ids
        

    #counts number of transmembrane parts in a protien
    def transmem_counter(self, s):
        if isinstance(s, str):
            return int(s.count("TRANSMEM"))
        else:
            return s

    #Isolates the n..n part of signal peptides
    #SIGNAL 1..22; /evidence="ECO:0000255 --> 1..22
    def isolateSignal(self, s):
        if not isinstance(s, float):
            num = re.search(r"\d+\.\.\d+", s)
            if num:
                return num.group()
            else:
                return None
        else:
            return None
        

    def downloadUniprot(self):
        #Checks for and downloads UniProtKB tsv file
        print("Looking for uniprot gene file")
        if os.path.exists(self.uniprot_file):
                print("TSV  File Found")

        else:
                print(f"Downloading {self.uniprot_file}")
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
        self.uniprot_genes = pd.read_csv(self.uniprot_file, sep='\t')
    
    def downloadDat(self):
        #Checks for and downloads UniProtKB .dat file
        print("Looking for uniprot dat file")
        gz_filename = self.isoformFile + ".gz"
        

        if os.path.exists(self.isoformFile):
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
                with open(self.isoformFile, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Decompressed: {self.isoformFile}")


    def uniprotMod(self):
        self.uniprot_genes['Protein existence'] = self.uniprot_genes['Protein existence'].apply(self.level_converter)
        self.uniprot_genes['Reviewed'] = self.uniprot_genes['Reviewed'].apply(self.reviewed)
        self.uniprot_genes['Transmembrane'] = self.uniprot_genes['Transmembrane'].apply(self.transmem_counter)
        self.uniprot_genes['Signal peptide'] = self.uniprot_genes['Signal peptide'].apply(self.isolateSignal)

    def datParser(self):
        with open(self.isoformFile) as UPfile:
            for line in UPfile:
                    if "AC   " in line:
                        changing_keys = line.strip().replace(" ","").replace("AC","").split(";")   
                        for key in changing_keys:
                            self.key_words[key] = [] 
                                
                    if "KW   " in line:
                        words = line.strip().replace(" ","").replace("KW","").replace(".",";")

                        for key in changing_keys:
                            if key in self.key_words:
                                self.key_words[key].append(words) 

                    if "Sequence=Displayed;" in line:
                            id = line.strip().replace(" ","").replace("CC","").split(";")[0].split("=")[1]
                            uniProtID = id.split("-")[0]
                            self.isoformDict[uniProtID] = id.split(",")[0]
                
                    if "RefSeq; NP_" in line:
                            isoNum = re.search(r"\[(.*?)\]", line)
                            if isoNum:
                                    isoNum = isoNum.group().replace("[","").replace("]","")
                                    prefix = isoNum.split("-")[0]
                                    if prefix in self.isoformDict and isoNum  == self.isoformDict[prefix] and prefix not in self.refSeqDict:
                                            self.refSeqDict[isoNum] = re.search(r"NP_(.*?);", line).group().replace(";","")
                    elif "RefSeq; XP_" in line:
                            isoNum = re.search(r"\[(.*?)\]", line)
                            if isoNum:
                                    isoNum = isoNum.group().replace("[","").replace("]","")
                                    prefix = isoNum.split("-")[0]
                                    if prefix in self.isoformDict and isoNum  == self.isoformDict[prefix] and prefix not in self.refSeqDict:
                                            self.refSeqDict[isoNum] = re.search(r"XP_(.*?);", line).group().replace(";","")
                    elif "RefSeq; NM_" in line:
                            isoNum = re.search(r"\[(.*?)\]", line)
                            if isoNum:
                                    isoNum = isoNum.group().replace("[","").replace("]","")
                                    prefix = isoNum.split("-")[0]
                                    if prefix in self.isoformDict and isoNum  == self.isoformDict[prefix] and prefix not in self.refSeqDict:
                                            self.refSeqDict[isoNum] = re.search(r"NM_(.*?);", line).group().replace(";","")




    def uniprotParse(self):
        #Creates a list of all the genes from GENCODE and combines them with their IDs, memory created to hold data collected from UniProtKB
        id_list = self.gene_df['gene_id'].tolist()
        name_list = self.gene_df['gene_name'].tolist()
        for i in range(len(id_list)):
            self.gene_dict[id_list[i]] = {"gene_id": id_list[i], "genecode_name": name_list[i], "gencode_symbol": None, "ENSP": [], "ENST":[], "uniprot_id": [], "reviewed": [], "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": [], "isoform":[], "EC Number":[], "Num Transmembrane Regions":[], "Signal Peptide":[], "refSeq Number":[], "Key Words":[]}
        print("Number of genes from GENCODE:", len(self.gene_dict))
        self.uniprotMod()
        print("Making connections with gene ids")
        # Has correct UniprotKB ID and the ENSG Number, hand chosen gene
        exceptions_dict = {"Q6UXT6":"ENSG00000228336"}

        #collects data from UniProtKB
        count = 0
        rows = 0
        #Links GENCODE genes to UniProtKB IDs through ENSG numbers 
        print("Selected UniProtIDs that have no listed Cannonical Isoform")
        for index, row in self.uniprot_genes.iterrows():
            row_dict = self.clean_string(row['Ensembl'], row['Entry'])
            if len(row_dict) >= 2:
                rows += 1
            for i in range(0, len(row_dict)):
                gene = row_dict[i]["gene_id"]		
                if gene in self.gene_dict and row['Reviewed']:
                                self.gene_dict[gene]['reviewed'].append(row['Reviewed'])
                                self.gene_dict[gene]['entry_name'].append(row['Entry Name'])
                                self.gene_dict[gene]['uniprot_id'].append(row['Entry'])
                                self.gene_dict[gene]['description'].append(row['Protein names']) 
                                self.gene_dict[gene]['protein length'].append(row['Length'])
                                self.gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                                self.gene_dict[gene]['found_with'].append("gene_id")
                                self.gene_dict[gene]['evidence'].append(row['Protein existence'])
                                self.gene_dict[gene]['entry_type'].append(row['Reviewed'])
                                self.gene_dict[gene]['ENSP'].append(row_dict[i]["ensp"])
                                self.gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                                self.gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])

                                self.gene_dict[gene]['EC Number'].append(row["EC number"])
                                self.gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])   
                                self.gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                                self.gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"])
                                self.gene_dict[gene]['Key Words'].append(self.key_words.get(row['Entry']))

                                count +=1
                elif gene in self.gene_dict and True not in self.gene_dict[gene]['reviewed']:
                                self.gene_dict[gene]['reviewed'].append(row['Reviewed'])
                                self.gene_dict[gene]['entry_name'].append(row['Entry Name'])
                                self.gene_dict[gene]['uniprot_id'].append(row['Entry'])
                                self.gene_dict[gene]['description'].append(row['Protein names'])
                                self.gene_dict[gene]['protein length'].append(row['Length'])
                                self.gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                                self.gene_dict[gene]['found_with'].append("gene_id")
                                self.gene_dict[gene]['evidence'].append(row['Protein existence'])
                                self.gene_dict[gene]['entry_type'].append(row['Reviewed'])
                                self.gene_dict[gene]['ENSP'].append(row_dict[i]["ensp"])
                                self.gene_dict[gene]['ENST'].append(row_dict[i]["trans_id"])
                                self.gene_dict[gene]['isoform'].append(row_dict[i]["isoform"])
                                self.gene_dict[gene]['EC Number'].append(row["EC number"])
                                self.gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])
                                self.gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                                self.gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"]) 
                                self.gene_dict[gene]['Key Words'].append(self.key_words.get(row['Entry']))

            if row['Entry'] in exceptions_dict:
                        gene = exceptions_dict[row['Entry']]
                        self.gene_dict[gene]['reviewed'].append(row['Reviewed'])
                        self.gene_dict[gene]['entry_name'].append(row['Entry Name'])
                        self.gene_dict[gene]['uniprot_id'].append(row['Entry'])
                        self.gene_dict[gene]['description'].append(row['Protein names'])
                        self.gene_dict[gene]['protein length'].append(row['Length'])
                        self.gene_dict[gene]['gene_symbol'].append(row['Gene Names'])
                        self.gene_dict[gene]['found_with'].append("gene_id")
                        self.gene_dict[gene]['evidence'].append(row['Protein existence'])
                        self.gene_dict[gene]['entry_type'].append(row['Reviewed'])
                        self.gene_dict[gene]['ENSP'].append(None)
                        self.gene_dict[gene]['ENST'].append(None) 
                        self.gene_dict[gene]['isoform'].append(self.isoformDict.get(row['Entry']))
                        self.gene_dict[gene]['EC Number'].append(row["EC number"])
                        self.gene_dict[gene]['Num Transmembrane Regions'].append(row["Transmembrane"])
                        self.gene_dict[gene]['Signal Peptide'].append(row["Signal peptide"])
                        self.gene_dict[gene]['refSeq Number'].append(self.refSeqDict.get(self.isoformDict.get(row['Entry']))) 
                        self.gene_dict[gene]['Key Words'].append(self.key_words.get(row['Entry']))

        print("Making connections with gene symbols")
        gene_symbols_dict = {}
        for i in self.gene_dict:
            gene_symbols_dict[self.gene_dict[i]["genecode_name"]] = i
        symbol_count = 0

        #Looks for GENCODE genes' link to UniProtKB through the gene symbols
        for index, row in self.uniprot_genes.iterrows():
            ids = str(row['Gene Names']).replace(';', '').split(' ')
            for name in ids:
                if name.strip() in gene_symbols_dict:
                                ensg = gene_symbols_dict[name.strip()]

                                if self.gene_dict[ensg]['found_with'] == [] or (True not in self.gene_dict[ensg]['reviewed'] and row['Reviewed']): 
                                            #print(gene_dict[i]['gencode_symbol'], ids)
                                            self.gene_dict[ensg]['entry_name'].append(row['Entry Name'])
                                            self.gene_dict[ensg]['uniprot_id'].append(row['Entry'])
                                            self.gene_dict[ensg]['description'].append(row['Protein names'])
                                            self.gene_dict[ensg]['protein length'].append(row['Length'])
                                            self.gene_dict[ensg]['gene_symbol'].append(name.strip())
                                            self.gene_dict[ensg]['EC Number'].append(row["EC number"])
                                            self.gene_dict[ensg]['Num Transmembrane Regions'].append(row["Transmembrane"])
                                            self.gene_dict[ensg]['found_with'].append("gene_name")
                                            
                                            self.gene_dict[ensg]['reviewed'].append(row['Reviewed'])
                                            self.gene_dict[ensg]['evidence'].append(row['Protein existence'])
                                            self.gene_dict[ensg]['entry_type'].append(row['Reviewed'])

                                            #Code is added for next clean_data.py
                                            self.gene_dict[ensg]['ENSP'].append(None)
                                            self.gene_dict[ensg]['ENST'].append(None)
                                            self.gene_dict[ensg]['isoform'].append(self.isoformDict.get(row['Entry Name']))
                                            self.gene_dict[ensg]['Signal Peptide'].append(row["Signal peptide"])
                                            self.gene_dict[ensg]['refSeq Number'].append(self.refSeqDict.get(self.isoformDict.get(row['Entry'])))
                                            self.gene_dict[ensg]['Key Words'].append(self.key_words.get(row['Entry']))
                                            symbol_count += 1 
        extra = 0
        not_found = 0
        #Checks for empty ENSG numbers
        for i in self.gene_dict:
            if not self.gene_dict[i]['entry_name']:
                not_found += 1
            if True not in self.gene_dict[i]["reviewed"]:
                extra += 1

        print("Num of reviewed", count)
        print("Num of non_reviewed", extra)
        #print("Num found with symbol", symbol_count)
        print("Num empty", not_found)
        print("Number of rows with 2 or more:", rows)                    

    def tableBuilder(self):
        print(f"Making Table: {self.output_file}")
        final_frame = pd.DataFrame(self.gene_dict).T

        merged_df = pd.merge(self.gene_df, final_frame, on='gene_id', how='inner')

        merged_df.to_excel(self.output_file)
        print("File made")  
        
    def run(self):
         self.downloadUniprot()
         self.downloadDat()
         self.datParser()
         self.uniprotParse()
         self.tableBuilder()

if __name__ == "__main__":
    processor = UniProtProcessor()
    processor.run()
 

