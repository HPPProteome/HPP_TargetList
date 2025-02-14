#Libraries
import os
import requests
import zipfile
import pandas as pd
import subprocess
import sys

class ProteinAtlasProcessor():
    def __init__(self):
        self.gene_file = "cleaned_table.xlsx"
        self.gene_data = pd.read_excel(self.gene_file)

        self.rna_file = "rna_tissue_consensus.tsv.zip"
        self.output_file = "updatedPE.xlsx"


    def downloadFile(self):
        print("Looking for RNA_expression file")
        if os.path.exists(self.rna_file):
            print("RNA expressions File Found")
        else:
            
            print("Downloading", self.rna_file)
            url = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
            output_zip_file = "rna_tissue_consensus.tsv.zip"

            print("Downloading", url)
            attempt = 0
            max_attempt = 3
            while attempt < max_attempt:
            
                try:
                    response = requests.get(url, stream=True, timeout=10)
                    with open(output_zip_file, 'wb') as f:
                        f.write(response.content)
                    print("Downloaded", output_zip_file)
                    break
                except requests.exceptions.Timeout:
                    attempt += 1
                    print(f"Request timed out, trying again {attempt}/{max_attempt}")
                except requests.exceptions.RequestException as e:
                    print("Error occured:", e)
                    break

            if attempt == max_attempt:
                print("File failed to download")
                sys.exit("Exiting Program")

    def tissueInformationAdder(self):
        #Reads data 
        rna_data = pd.read_csv('rna_tissue_consensus.tsv.zip',compression='zip', sep='\t')
        

        #How high a nTPM score must be to assign a PE2
        threshold = 1 #change to be lower or higher

        self.gene_data['Suggested PE'] = ''
        self.gene_data['Highest nTPM Score'] = ''
        self.gene_data[f"Tissues with nTPM Score Above {threshold} (/50)"] = ''
        rna_dict = {}

        #Collects all nTPM scores for every gene present  
        for index, row in rna_data.iterrows():
            gene = row['Gene']
            if gene not in rna_dict:
                rna_dict[gene] = []
            rna_dict[gene].append(row['nTPM'])

        count = 0
        changed = 0

        #Checks to see if a PE > 2 meets the qualifications to be assigned 2.
        for index, row in self.gene_data.iterrows():
            if row['gene_id'] in rna_dict:
                count += 1
                
                self.gene_data.at[index, 'Highest nTPM Score'] = max(rna_dict[row['gene_id']])
                self.gene_data.at[index, f"Tissues with nTPM Score Above {threshold} (/50)"]  = len([score for score in rna_dict[row['gene_id']] if score >= threshold])		

                if max(rna_dict[row['gene_id']]) >= threshold and row['evidence'] > 2:
                    self.gene_data.at[index, 'Suggested PE'] = 2
                    changed += 1

        print("Number of genes in RNA tsv file:", count)
        print("Number of genes with new suggested PE score:", changed)

    def tableBuilder(self):

        print(f"Making Frame: {self.output_file}")
        self.gene_data.to_excel(self.output_file, index=False)
        print("done")

    def run(self):
        self.downloadFile()
        self.tissueInformationAdder()
        self.tableBuilder()

if __name__ == "__main__":
    processor = ProteinAtlasProcessor()
    processor.run()
