#Libraries
import pandas as pd
import os
import subprocess


class CleanDataProcessor:
    def __init__(self, output_file="cleaned_data.xlsx"):
        self.entries_file = "uniprot_output.xlsx"
        self.gene_data = pd.read_excel("uniprot_output.xlsx")
        self.columns = ['uniprot_id', 'entry_name', 'gene_symbol', 'description', 'protein length', 'entry_type', 'evidence', 'found_with', 'isoform', 'EC Number', 'Num Transmembrane Regions', 'refSeq Number', 'ENSP', 'ENST', 'Signal Peptide']

        self.strange_genes = []

    def makeList(self, collection):
        return collection.replace('[', '').replace(']', '').replace("'", '').strip().split(',')
    
    def stringListToList(self, input_string):
        return eval(input_string)




    def dataModifier(self):
        self.gene_data.columns = self.gene_data.columns.str.strip()

        #Since reading an excel file into a dataframe creates a single string, values  need to be re-seperated into lists
        for i in self.columns:
            self.gene_data[i] = self.gene_data[i].apply(self.makeList)

        #Needs to be specially made into list
        self.gene_data["Key Words"] = self.gene_data["Key Words"].apply(self.stringListToList)


    def falseOverTrue(self):
        #Optimizes for using gene_id over gene_symbol. Catches any True gene_symbols and False gene_id enteries to store in diff file.
        #Will keep False gene_id over True gene_symbol entry.

        self.gene_data['Reviewed Entry Available'] = ''
        for index, row in self.gene_data.iterrows():
            found_type = row['found_with']
            contains_id = 'gene_id' in [val.strip() for val in found_type]
            

            if contains_id:
                found_right = [types.strip() == 'gene_id' for types in found_type]

                self.gene_data.at[index, 'uniprot_id'] = [row['uniprot_id'][i] for i in range(len(found_right)) if found_right[i]]  
                self.gene_data.at[index, 'entry_name'] = [row['entry_name'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'gene_symbol'] = [row['gene_symbol'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'description'] = [row['description'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'protein length'] = [row['protein length'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'evidence'] = [row['evidence'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'entry_type'] = [row['entry_type'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'ENSP'] = [row['ENSP'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'ENST'] = [row['ENST'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'isoform'] = [row['isoform'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'Num Transmembrane Regions'] = [row['Num Transmembrane Regions'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'EC Number'] = [row['EC Number'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'found_with'] = [row['found_with'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'refSeq Number'] = [row['refSeq Number'][i] for i in range(len(found_right)) if found_right[i]]
                self.gene_data.at[index, 'Key Words'] = [row['Key Words'][i] for i in range(len(found_right)) if found_right[i]]
                for i in range(len(found_type)):
                    if row['found_with'][i].strip()  == "gene_name" and row['entry_type'][i].strip():
                    
                        self.strange_genes.append({col: row[col][i] for col in row.index if isinstance(row[col], list)})
                        self.gene_data.at[index, 'Reviewed Entry Available'] = 'yes'

        
        self.strange_genes = pd.DataFrame(self.strange_genes, columns=self.columns).fillna("")
        

    def cleanData(self):
        #Gets rid of any false entries if a GENCODE gene has a corosponding reviewed UniProt entry
        for index, row in self.gene_data.iterrows():
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

                self.gene_data.at[index, 'uniprot_id'] = filtered_uniprot_ids
                self.gene_data.at[index, 'entry_name'] = filtered_entry_names
                self.gene_data.at[index, 'gene_symbol'] = filtered_gene_symbols
                self.gene_data.at[index, 'description'] = filtered_descriptions
                self.gene_data.at[index, 'protein length'] = filtered_protein_lengths
                self.gene_data.at[index, 'evidence'] = filtered_evidence
                self.gene_data.at[index, 'entry_type'] = filtered_types
                
        
                self.gene_data.at[index, 'ENSP'] = [row['ENSP'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'ENST'] = [row['ENST'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'isoform'] = [row['isoform'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]

                self.gene_data.at[index, 'Num Transmembrane Regions'] = [row['Num Transmembrane Regions'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'EC Number'] = [row['EC Number'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'found_with'] = [row['found_with'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'refSeq Number'] = [row['refSeq Number'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
                self.gene_data.at[index, 'Key Words'] = [row['Key Words'][i] for i in range(len(entry_types_bool)) if entry_types_bool[i]]
        #Makes sure that only the lowest  level of exsistance proteins are kept [1,1,4] --> [1,1]
        for index, row in self.gene_data.iterrows():
            if len(row['evidence']) > 1:
                level = [int(exist) for exist in row['evidence']]
                best = min(level)
                keeper = []
                for i in level:
                    keeper.append(i==best)
                self.gene_data.at[index, 'uniprot_id'] = [row['uniprot_id'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'entry_name'] = [row['entry_name'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'gene_symbol'] = [row['gene_symbol'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'description'] = [row['description'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'protein length'] = [row['protein length'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'evidence'] = [row['evidence'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'entry_type'] = [row['entry_type'][i].strip() for i in range(len(keeper)) if keeper[i]]
                
                
                self.gene_data.at[index, 'ENSP'] = [row['ENSP'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'ENST'] = [row['ENST'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'isoform'] = [row['isoform'][i].strip() for i in range(len(keeper)) if keeper[i]]

                self.gene_data.at[index, 'EC Number'] = [row['EC Number'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'found_with'] = [row['found_with'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'Signal Peptide'] = [row['Signal Peptide'][i].strip() for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'refSeq Number'] = [row['refSeq Number'][i] for i in range(len(keeper)) if keeper[i]]
                self.gene_data.at[index, 'Key Words'] = [row['Key Words'][i] for i in range(len(keeper)) if keeper[i]]



        #Uses the entry with the connocial ENSP and ENST numbers (if availible), by taking the only avaible ENSP number
        for index, row in self.gene_data.iterrows():
            if isinstance(row['ENSP'], list) and len(row['ENSP']) > 1:
                rowSet = set(row['ENSP'])
                rowSet.discard('None')
                if len(rowSet) == 1:
                        ENSP = list(rowSet)[0]
                        i = row['ENSP'].index(ENSP)
                        
                        self.gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
                        self.gene_data.at[index, 'entry_name'] = row['entry_name'][i]
                        self.gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]
                        self.gene_data.at[index, 'description'] = row['description'][i]
                        self.gene_data.at[index, 'protein length'] = row['protein length'][i]
                        self.gene_data.at[index, 'evidence'] = row['evidence'][i]
                        self.gene_data.at[index, 'entry_type'] = row['entry_type'][i]
                        self.gene_data.at[index, 'ENSP'] = row['ENSP'][i]
                        self.gene_data.at[index, 'ENST'] = row['ENST'][i]
                        self.gene_data.at[index, 'isoform'] = row['isoform'][i]
                        self.gene_data.at[index, 'EC Number'] = row['EC Number'][i]
                        self.gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][i]
                        self.gene_data.at[index, 'found_with'] = row['found_with'][i]
                        self.gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][i]
                        self.gene_data.at[index, 'refSeq Number'] = row['refSeq Number'][i]
                        self.gene_data.at[index, 'Key Words'] = row['Key Words'][i]
        





        repeats = 0
        #Uses 1 entry when there are muiltple replicas
        for index, row in self.gene_data.iterrows():
                if isinstance(row['uniprot_id'], list) and len(row['uniprot_id']) > 1:
                    stripped_row = [length.strip() for length in row['uniprot_id']]
                    if len(set(stripped_row)) == 1:
                            repeats += 1
                            self.gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][0]
                            self.gene_data.at[index, 'entry_name'] = row['entry_name'][0]
                            self.gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][0]
                            self.gene_data.at[index, 'description'] = row['description'][0]
                            self.gene_data.at[index, 'protein length'] = row['protein length'][0]
                            self.gene_data.at[index, 'evidence'] = row['evidence'][0]
                            self.gene_data.at[index, 'entry_type'] = row['entry_type'][0]

                            self.gene_data.at[index, 'ENSP'] = row['ENSP'][0]
                            self.gene_data.at[index, 'ENST'] = row['ENST'][0]
                            self.gene_data.at[index, 'isoform'] = row['isoform'][0]
                            self.gene_data.at[index, 'EC Number'] = row['EC Number'][0]
                            self.gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][0]
                            self.gene_data.at[index, 'found_with'] = row['found_with'][0]
                            self.gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][0]
                            self.gene_data.at[index, 'refSeq Number'] = row['refSeq Number'][0]
                            self.gene_data.at[index, 'Key Words'] = row['Key Words'][0]
        print("Number of repeats:", repeats)

        #Chooses highest CDS length when muiltiple entries could be a possible match
        for index, row in self.gene_data.iterrows():
            if len(row['entry_type']) > 1 and isinstance(row['entry_type'], list):
                mini_dict = {"longest":0, "index":None}
                for i in range(0, len(row['protein length'])):
                    if int(row['protein length'][i]) > mini_dict["longest"]:
                        mini_dict["longest"] = int(row['protein length'][i])
                        mini_dict["index"] = i
            
                i = mini_dict["index"]
                self.gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][i]
                self.gene_data.at[index, 'entry_name'] = row['entry_name'][i]
                self.gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][i]
                self.gene_data.at[index, 'description'] = row['description'][i]
                self.gene_data.at[index, 'protein length'] = row['protein length'][i]
                self.gene_data.at[index, 'evidence'] = row['evidence'][i]
                self.gene_data.at[index, 'entry_type'] = row['entry_type'][i]
                self.gene_data.at[index, 'ENSP'] = row['ENSP'][i]
                self.gene_data.at[index, 'ENST'] = row['ENST'][i]
                self.gene_data.at[index, 'isoform'] = row['isoform'][i]
                self.gene_data.at[index, 'EC Number'] = row['EC Number'][i]
                self.gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][i]
                self.gene_data.at[index, 'found_with'] = row['found_with'][i]
                self.gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][i]
                self.gene_data.at[index, 'refSeq Number'] = row['refSeq Number'][i]
                self.gene_data.at[index, 'Key Words'] = row['Key Words'][i]

        
        #Cleans data so that is no longer a list
        for index, row in self.gene_data.iterrows():
            if isinstance(row['entry_name'], list) and len(row['entry_name']) == 1:
                            self.gene_data.at[index, 'uniprot_id'] = row['uniprot_id'][0]
                            self.gene_data.at[index, 'entry_name'] = row['entry_name'][0]
                            self.gene_data.at[index, 'gene_symbol'] = row['gene_symbol'][0]
                            self.gene_data.at[index, 'description'] = row['description'][0]
                            self.gene_data.at[index, 'protein length'] = row['protein length'][0]
                            self.gene_data.at[index, 'evidence'] = row['evidence'][0]
                            self.gene_data.at[index, 'entry_type'] = row['entry_type'][0]
                            self.gene_data.at[index, 'ENSP'] = row['ENSP'][0].strip()
                            self.gene_data.at[index, 'ENST'] = row['ENST'][0].strip()
                            self.gene_data.at[index, 'EC Number'] = row['EC Number'][0]
                            self.gene_data.at[index, 'Num Transmembrane Regions'] = row['Num Transmembrane Regions'][0]
                            self.gene_data.at[index, 'found_with'] = row['found_with'][0]
                            self.gene_data.at[index, 'Signal Peptide'] = row['Signal Peptide'][0]
                            self.gene_data.at[index, 'refSeq Number'] = row['refSeq Number'][0]
                            #self.gene_data.at[index, 'Key Words'] = row['Key Words'][0]
                            if row['isoform'][0] != None and not row['isoform'][0].startswith("ENSG"):
                                self.gene_data.at[index, 'isoform'] = row['isoform'][0].strip()
                            else:
                                self.gene_data.at[index, 'isoform'] = 'None'


                            if row['Key Words'] and None not in row['Key Words']:
                                if isinstance(row['Key Words'][0], list):
                                    self.gene_data.at[index, 'Key Words'] = ", ".join(row['Key Words'][0]).replace(";", ", ")
                                else:
                                    self.gene_data.at[index, 'Key Words'] = ", ".join(row['Key Words']).replace(";", ", ")
                            else:
                                self.gene_data.at[index, 'Key Words'] = ""    
                            self.gene_data.at[index, 'Key Words'] = str(self.gene_data.at[index, 'Key Words']).replace(";", " ,").replace("[", "").replace("]", "").replace(", ,", ",")          


    def makeTable(self):

        print(self.gene_data.shape)

        print("Making Frame: cleaned_table.xlsx and False_ensg_over_name.xlsx")
        self.gene_data.to_excel("cleaned_table.xlsx", index=False)
        self.strange_genes.to_excel("False_ensg_over_name.xlsx", index=False)
        print("done")



    def run(self):
        self.dataModifier()
        self.falseOverTrue()
        self.cleanData()
        self.makeTable()




if __name__ == "__main__":
    processor = CleanDataProcessor()
    processor.run()    
    
