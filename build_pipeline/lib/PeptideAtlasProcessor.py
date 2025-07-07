import requests
import os
import pandas as pd
import subprocess
import sys

class PeptideAtlasProcessor:
    def __init__(self):
        self.PeptideAtlas_core_proteome_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetCoreProteomeMapping?mapping_id=116&atlas_build_id=592&redundancy_constraint=on&QUERY_NAME=AT_GetCoreProteomeMapping&apply_action=QUERY&output_mode=tsv"
        self.PeptideAtlas_all_proteins_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins?row_limit=5000000&atlas_build_id=592&organism_id=2&apply_action=QUERY&output_mode=tsv"
        self.core_atlas_file  = "peptideAtlas.tsv"
        self.all_atlas_file  = "Secoundary_peptideAtlas.tsv"
        self.gene_file = "updatedPE.xlsx"
        self.gene_df = None
        self.atlas_dict = {}
        self.secound_atlas_dict = {}
        
    def download_file(self, url, filename):
        attempt, max_attempts = 0, 3
        while attempt < max_attempts:
            try:
                print(f"Downloading {filename} from Peptide Atlas")
                response = requests.get(url, stream=True, timeout=30)
                response.raise_for_status()
                with open(filename, "wb") as file:
                    file.write(response.content)
                print(f"File downloaded and saved as {filename}")
                return
            except requests.exceptions.RequestException as e:
                print(f"An error occurred: {e}")
            except requests.exceptions.Timeout:
                attempt += 1
                print(f"Request Timed Out, trying again {attempt}/{max_attempts}")
        print(f"Failed to download {filename}")
        sys.exit("Exiting program")
        
    def ensure_files_exist(self):
        if not os.path.exists(self.core_atlas_file):
            self.download_file(self.PeptideAtlas_core_proteome_url, self.core_atlas_file)
        else:
            print("Peptide Atlas file found")
        
        if not os.path.exists(self.all_atlas_file):
            self.download_file(self.PeptideAtlas_all_proteins_url, self.all_atlas_file)
        else:
            print("Secondary Atlas file found")
        
        print(f"{self.gene_file} found")
        
    def load_data(self):
        self.gene_df = pd.read_excel(self.gene_file)
        atlas_df = pd.read_csv(self.core_atlas_file, sep='\t')
        secound_atlas_df = pd.read_csv(self.all_atlas_file, sep='\t', low_memory=False)
        
        for _, row in atlas_df.iterrows():
            if isinstance(row['Ensembl_Accession'], str):
                self.atlas_dict[row['Ensembl_Accession']] = {
                    "catagory": row['PeptideAtlas_Category'],
                    "observed": row['nobs'],
                    "distinct": row['npep'],
                    "unique": row['nunipep'],
                    "uniprot": row['accession']
                }
        
        for _, row in secound_atlas_df.iterrows():
            if isinstance(row['biosequence_name'], str):
                self.secound_atlas_dict[row['biosequence_name']] = {
                    "catagory": row['presence_level'],
                    "observed": row['n_observations'],
                    "distinct": row['n_distinct_peptides']
                }
        
    def process_data(self):
        self.gene_df['PeptideAtlas Category'] = ''
        self.gene_df['Observed'] = ''
        self.gene_df['Distinct'] = ''
        self.gene_df['Uniquely Mapping'] = ''
        
        not_in, mismatch, secound_dict = 0, 0, 0
        
        for index, row in self.gene_df.iterrows():
            gene, uniprot_id = row['gene_id'], row['uniprot_id']
            if gene in self.atlas_dict:
                self.gene_df.at[index, 'PeptideAtlas Category'] = self.atlas_dict[gene]['catagory']
                self.gene_df.at[index, 'Observed'] = self.atlas_dict[gene]['observed']
                self.gene_df.at[index, 'Distinct'] = self.atlas_dict[gene]['distinct']
                self.gene_df.at[index, 'Uniquely Mapping'] = self.atlas_dict[gene]['unique']
                if self.atlas_dict[gene]['uniprot'] != uniprot_id:
                    mismatch += 1
            elif uniprot_id in self.secound_atlas_dict:
                secound_dict += 1
                self.gene_df.at[index, 'PeptideAtlas Category'] = self.secound_atlas_dict[uniprot_id]['catagory']
                self.gene_df.at[index, 'Observed'] = self.secound_atlas_dict[uniprot_id]['observed']
                self.gene_df.at[index, 'Distinct'] = self.secound_atlas_dict[uniprot_id]['distinct']
            else:
                not_in += 1
        
        print(f"There are {mismatch} mismatched uniprot IDs")
        print(f"There are {secound_dict} genes found in the secondary atlas file")
        print(f"There are {not_in} genes not in Peptide Atlas")
        
        self.gene_df.to_excel("atlasLink.xlsx", index=False)
        print("Processing complete. File saved as atlasLink.xlsx")
        
    def run(self):
        self.ensure_files_exist()
        self.load_data()
        self.process_data()

if __name__ == "__main__":
    processor = PeptideAtlasProcessor()
    processor.run()

