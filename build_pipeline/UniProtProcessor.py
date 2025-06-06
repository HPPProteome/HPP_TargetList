import os
import requests
import shutil
import gzip
import pandas as pd
import subprocess
import re
import sys
from Bio import SeqIO


class UniProtProcessor:
    def __init__(self):
        self.uniprot_file = "uniprot.tsv.gz"
        self.gene_file = "protien_coding_genes.xlsx"
        self.isoformFile = "Uniprot.dat"
        self.additional_isoformFile = "additional_datFile.dat"
        
        self.output_file = "uniprot_output.xlsx"
        
        self.gene_df = pd.read_excel(self.gene_file)
        self.uniprot_genes = pd.read_csv(self.uniprot_file, sep='\t')
        
        self.gene_dict = {}
        self.key_words = {} #UniProtID (s) : Associated tags
        self.isoformDict = {} # UniProtID : Cannonical Isoform
        self.refSeqDict = {}  #Isoform : refSequence


        self.gencode_fasta = "gencode.pc_translations.fa.gz"
        self.uniprot_fasta = "uniprot.fa"
        self.added = 0
        self.gencode_seq = {} #Stored ensgNum{sequence:ENSP}
        self.uniprot_seq = {} #Stores sequence:ID

        self.exceptions_file = "Manual_ENSG_UP_associations.tsv"
        self.exceptions_df = pd.read_csv('Manual_ENSG_UP_associations.tsv', sep='\t')
        self.exceptions_dict = {}

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
        
    #Creates a dictionary with the ensamble row, then picks the isoform(if present) found in Uniprot.dat and gets the ensp, enst and ensg number
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
        
    def readExceptionFile(self):
        
        for index, row in self.exceptions_df.iterrows():
            if row["UP id"] not in self.exceptions_dict:
                self.exceptions_dict[row["UP id"]] = []
            self.exceptions_dict[row["UP id"]].append(row["ENSG identifier"])


    def uniprotMod(self):
        self.uniprot_genes['Protein existence'] = self.uniprot_genes['Protein existence'].apply(self.level_converter)
        self.uniprot_genes['Reviewed'] = self.uniprot_genes['Reviewed'].apply(self.reviewed)
        self.uniprot_genes['Transmembrane'] = self.uniprot_genes['Transmembrane'].apply(self.transmem_counter)
        self.uniprot_genes['Signal peptide'] = self.uniprot_genes['Signal peptide'].apply(self.isolateSignal)

    def datParser_function(self, file):
        with open(file) as UPfile:
            for line in UPfile:
                    if "ID   " in line:
                        changing_keys = []

                    if "AC   " in line:
                        accession_data = line[5:].strip()  # skips "AC   "
                        key_list = [k.strip() for k in accession_data.split(";") if k.strip()]                        

                        for key in key_list:
                            if key not in self.key_words:
                                self.key_words[key] = {"Key Words": [], "PPI":0, "Main ID":""} #Creates a space to hold info about key words and number of protien protien interaction 
                            if key not in changing_keys:
                                changing_keys.append(key)
                        

                    elif "KW  " in line:
                        words = line.strip().replace(" ","").replace("KW","").replace(".",";")

                        for key in changing_keys:
                            if key in self.key_words:
                                self.key_words[key]["Key Words"].append(words)
                            
                            if self.key_words[key]["Main ID"] == "":
                                self.key_words[key]["Main ID"] = changing_keys[0] 

                    if "IntAct=" in line:
                        for key in changing_keys:
                            self.key_words[key]["PPI"] += 1


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
    def datParser(self):
        self.datParser_function(self.isoformFile)
        self.datParser_function(self.additional_isoformFile)

    def sequenceLoader(self):
         
    
        for record in SeqIO.parse(self.uniprot_fasta, "fasta"):
                    header_parts = record.description.split('|')
                    id = header_parts[1]
                    seq = str(record.seq)
                    self.uniprot_seq[seq] = id

        with gzip.open(self.gencode_fasta, "rt") as unzipped_fasta:
                    for record in SeqIO.parse(unzipped_fasta, "fasta"):
                        header_parts = record.description.split('|')
                        ensp_number = header_parts[0]
                        ensg_number = header_parts[2].split(".")[0]
                        enst_number = header_parts[1]
                        seq = str(record.seq)

                        if ensg_number not in self.gencode_seq:
                            self.gencode_seq[ensg_number] = {}
                        self.gencode_seq[ensg_number][seq] = [ensp_number, enst_number]
        print("Files read")
        print("\nENSG numbers assigned:")


    def uniprotParse(self):
        self.readExceptionFile()
        #Creates a list of all the genes from GENCODE and combines them with their IDs, memory created to hold data collected from UniProtKB
        id_list = self.gene_df['gene_id'].tolist()
        name_list = self.gene_df['gene_name'].tolist()
        for i in range(len(id_list)):
            self.gene_dict[id_list[i]] = {"gene_id": id_list[i], "genecode_name": name_list[i], "gencode_symbol": None, "ENSP": [], "ENST":[], "uniprot_id": [], "reviewed": [], "entry_name": [], "gene_symbol": [], "description": [], "protein length": [], "entry_type": [], "evidence": [], "found_with": [], "isoform":[], "EC Number":[], "Num Transmembrane Regions":[], "Signal Peptide":[], "refSeq Number":[], "Key Words":[], "PPI":[], "Main ID":[]}
        print("Number of genes from GENCODE:", len(self.gene_dict))
        self.uniprotMod()
        print("Making connections with gene ids")
        print(f"Exceptions dict should be used {len(self.exceptions_dict)} times")
        # Has correct UniprotKB ID and the ENSG Number, hand chosen gene

        #collects data from UniProtKB
        count = 0
        rows = 0
        #Links GENCODE genes to UniProtKB IDs through ENSG numbers 
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

                                self.gene_dict[gene]['EC Number'].append(row.get("EC number", None))
                                self.gene_dict[gene]['Num Transmembrane Regions'].append(row.get("Transmembrane", "None"))   
                                self.gene_dict[gene]['Signal Peptide'].append(row.get("Signal peptide", "None"))
                                self.gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"])
                                self.gene_dict[gene]['Key Words'].append(self.key_words.get(row['Entry'], {}).get("Key Words", "Not Present"))
                                self.gene_dict[gene]['PPI'].append(self.key_words.get(row['Entry'], {}).get("PPI", "Not Present"))
                                self.gene_dict[gene]['Main ID'].append(self.key_words.get(row['Entry'], {}).get("Main ID", "Not Present"))


                                count +=1
                elif gene in self.gene_dict and (True not in self.gene_dict[gene]['reviewed'] or ('gene_id' in self.gene_dict[gene]['found_with'] and row['Reviewed'])):
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
                                self.gene_dict[gene]['EC Number'].append(row.get("EC number", None))
                                self.gene_dict[gene]['Num Transmembrane Regions'].append(row.get("Transmembrane", "None"))
                                self.gene_dict[gene]['Signal Peptide'].append(row.get("Signal peptide", "None"))
                                self.gene_dict[gene]['refSeq Number'].append(row_dict[i]["refSeq"]) 
                                self.gene_dict[gene]['Key Words'].append(self.key_words.get(row['Entry'], {}).get("Key Words", "Not Present"))
                                self.gene_dict[gene]['PPI'].append(self.key_words.get(row['Entry'], {}).get("PPI", "Not Present"))
                                self.gene_dict[gene]['Main ID'].append(self.key_words.get(row['Entry'], {}).get("Main ID", "Not Present"))

            if row['Entry'] in self.exceptions_dict:
                        print("exceptions_dict used")
                        genes = self.exceptions_dict[row['Entry']]
                        for gene in genes:
                            if gene in self.gene_dict:
                                self.gene_dict[gene]['reviewed'] = [row['Reviewed']]
                                self.gene_dict[gene]['entry_name'] = [row['Entry Name']]
                                self.gene_dict[gene]['uniprot_id'] = [row['Entry']]
                                self.gene_dict[gene]['description'] = [row['Protein names']]
                                self.gene_dict[gene]['protein length'] = [row['Length']]
                                self.gene_dict[gene]['gene_symbol'] = [row['Gene Names']]
                                self.gene_dict[gene]['found_with'] = ["Hand Selected"]
                                self.gene_dict[gene]['evidence'] = [row['Protein existence']]
                                self.gene_dict[gene]['entry_type'] = [row['Reviewed']]
                                self.gene_dict[gene]['ENSP'] = [None]
                                self.gene_dict[gene]['ENST'] = [None] 
                                self.gene_dict[gene]['isoform'] = [self.isoformDict.get(row['Entry'])]
                                self.gene_dict[gene]['EC Number'] = [row.get("EC number", None)]
                                self.gene_dict[gene]['Num Transmembrane Regions'] = [row.get("Transmembrane", "None")]
                                self.gene_dict[gene]['Signal Peptide'] = [row.get("Signal peptide", "None")]
                                self.gene_dict[gene]['refSeq Number'] = [self.refSeqDict.get(self.isoformDict.get(row['Entry']))]
 
                                self.gene_dict[gene]['Key Words'] = [self.key_words.get(row['Entry'], {}).get("Key Words", "Not Present")]
                                self.gene_dict[gene]['PPI'] = [self.key_words.get(row['Entry'], {}).get("PPI", "Not Present")]
                                self.gene_dict[gene]['Main ID'] = [self.key_words.get(row['Entry'], {}).get("Main ID", "Not Present")]
                            else:
                                print(f"There is no need for an exception for {gene}. It is not in the release")
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

                                if self.gene_dict[ensg]['found_with'] == [] or (True not in self.gene_dict[ensg]['reviewed'] and row['Reviewed']) or (row['Reviewed'] and "gene_name" in self.gene_dict[ensg]['found_with']): 
                                            #print(gene_dict[i]['gencode_symbol'], ids)
                                            self.gene_dict[ensg]['entry_name'].append(row['Entry Name'])
                                            self.gene_dict[ensg]['uniprot_id'].append(row['Entry'])
                                            self.gene_dict[ensg]['description'].append(row['Protein names'])
                                            self.gene_dict[ensg]['protein length'].append(row['Length'])
                                            self.gene_dict[ensg]['gene_symbol'].append(name.strip())
                                            self.gene_dict[ensg]['EC Number'].append(row.get("EC number", "None"))
                                            self.gene_dict[ensg]['Num Transmembrane Regions'].append(row.get("Transmembrane", "None"))
                                            self.gene_dict[ensg]['found_with'].append("gene_name")
                                            
                                            self.gene_dict[ensg]['reviewed'].append(row['Reviewed'])
                                            self.gene_dict[ensg]['evidence'].append(row['Protein existence'])
                                            self.gene_dict[ensg]['entry_type'].append(row['Reviewed'])

                                            #Code is added for next clean_data.py
                                            self.gene_dict[ensg]['ENSP'].append(None)
                                            self.gene_dict[ensg]['ENST'].append(None)
                                            self.gene_dict[ensg]['isoform'].append(self.isoformDict.get(row['Entry Name']))
                                            self.gene_dict[ensg]['Signal Peptide'].append(row.get("Signal peptide", None))
                                            self.gene_dict[ensg]['refSeq Number'].append(self.refSeqDict.get(self.isoformDict.get(row['Entry'])))
                                            

                                            self.gene_dict[ensg]['Key Words'].append(self.key_words.get(row['Entry'], {}).get("Key Words", "Not Present"))
                                            self.gene_dict[ensg]['PPI'].append(self.key_words.get(row['Entry'], {}).get("PPI", "Not Present"))
                                            self.gene_dict[ensg]['Main ID'].append(self.key_words.get(row['Entry'], {}).get("Main ID", "Not Present"))
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
        self.gene_table = pd.DataFrame(self.gene_dict).T
        self.merged_df = pd.merge(self.gene_df, self.gene_table, on='gene_id', how='inner')

    def SequenceLink(self):
        print("\nMaking connections via amino acid sequences\n")
        self.sequenceLoader()
        for index, row in self.merged_df.iterrows():
            if not row['uniprot_id'] and row['gene_id'] in self.gencode_seq:
                possible_seq = self.gencode_seq[row['gene_id']]
                for sequence in possible_seq:
                    if sequence in self.uniprot_seq:
                        uniprot_id = self.uniprot_seq[sequence]
                        ensp = possible_seq[sequence][0]
                        enst = possible_seq[sequence][1]
                        self.merged_df.at[index, 'uniprot_id'] = [uniprot_id]
                        self.merged_df.at[index, 'ENSP'] = [ensp]
                        self.merged_df.at[index, 'ENST'] = [enst]
                        
                        information = self.uniprot_genes[self.uniprot_genes['Entry'] == uniprot_id]

                        self.merged_df.at[index, 'reviewed'] = [information['Reviewed'].values[0]]
                        self.merged_df.at[index, 'entry_name'] = [information['Entry Name'].values[0]]
                        self.merged_df.at[index, 'gene_symbol'] = [information['Gene Names'].values[0]]
                        self.merged_df.at[index, 'description'] = [information['Protein names'].values[0]]
                        self.merged_df.at[index, 'protein length'] = [information['Length'].values[0]]

                        self.merged_df.at[index, 'entry_type'] = [information['Reviewed'].values[0]]
                        self.merged_df.at[index, 'evidence'] = [information['Protein existence'].values[0]]

                        self.merged_df.at[index, 'found_with'] = ["Sequence"]
                        self.merged_df.at[index, 'isoform'] = [self.isoformDict.get(uniprot_id)]
            

                        self.merged_df.at[index, 'EC Number'] = [information["EC number"].values[0]]
                        self.merged_df.at[index, 'Num Transmembrane Regions'] = [information["Transmembrane"].values[0]]
                        self.merged_df.at[index, 'Signal Peptide'] = [information["Signal peptide"].values[0]]
                        self.merged_df.at[index, 'refSeq Number'] = [self.refSeqDict.get(self.isoformDict.get(uniprot_id))]



                        self.merged_df.at[index, 'Key Words'] = [self.key_words.get(uniprot_id, {}).get("Key Words", "Not Present")]
                        self.merged_df.at[index, 'PPI'] = [self.key_words.get(uniprot_id, {}).get("PPI", "Not Present")]
                        self.merged_df.at[index, 'Main ID'] = [self.key_words.get(uniprot_id, {}).get("Main ID", "Not Present")]                              

                        self.added += 1
                        print(f"{row['gene_id']}: {uniprot_id}")
                        
                        break

    def tableBuilder(self):
        print(f"Making Table: {self.output_file}")

        

        self.merged_df.to_excel(self.output_file)
        print("File made")  
        
    def run(self):
         print("Parsing Files")
         self.datParser()
         print("Matching ENSG with UniprotKB IDs through Ensamble")
         self.uniprotParse()
         self.SequenceLink()
         self.tableBuilder()


if __name__ == "__main__":
    processor = UniProtProcessor()
    processor.run()
 

