import requests
import os
import pandas as pd
import subprocess
import sys

class PeptideAtlasProcessor:
    def __init__(self):
        self.atlas_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetCoreProteomeMapping/query_zsun_20241030-091733-695.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_zsun_20241030-091733-695&rs_page_size=1000000&output_mode=tsv"
        self.backup_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins/query_edeutsch_20250117-134052-808.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_edeutsch_20250117-134052-808&rs_page_size=1000000&output_mode=tsv"
        self.atlas_file = "peptideAtlas.tsv"
        self.backup_atlas_file = "Secoundary_peptideAtlas.tsv"
        self.gene_file = "updatedPE.xlsx"
        self.gene_df = None
        self.atlas_dict = {}
        self.secound_atlas_dict = {}
        
    def download_file(self, url, filename):
        attempt, max_attempts = 0, 3
        while attempt < max_attempts:
            try:
                print(f"Downloading {filename} from Peptide Atlas")
                response = requests.get(url, stream=True, timeout=10)
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
        if not os.path.exists(self.atlas_file):
            self.download_file(self.atlas_url, self.atlas_file)
        else:
            print("Peptide Atlas file found")
        
        if not os.path.exists(self.backup_atlas_file):
            self.download_file(self.backup_url, self.backup_atlas_file)
        else:
            print("Secondary Atlas file found")
        
        if not os.path.exists(self.gene_file):
            print(f"Can't find {self.gene_file}, running update_PE.py")
            try:
                subprocess.run(['python3', 'update_PE.py'], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running update_PE.py: {e}")
        print(f"{self.gene_file} found")
        
    def load_data(self):
        self.gene_df = pd.read_excel(self.gene_file)
        atlas_df = pd.read_csv(self.atlas_file, sep='\t')
        secound_atlas_df = pd.read_csv(self.backup_atlas_file, sep='\t', low_memory=False)
        
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

