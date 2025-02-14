#Importing Libraries

from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import sys

class GENCODEProcessor:
    def __init__(self, version):
        self.version = version
        self.gtf_file = "gencode.annotation.gtf.gz"
        self.output_file = "protien_coding_genes.xlsx"
        self.unexcepted_tags = ['readthrough_transcript', 'EnsEMBL_merge_exception', 'readthrough_gene', 'exp_conf', 'RNA_Seq_supported_only']
        self.unexcepted_trans_type = ['nonsense_mediated_decay', 'protein_coding_CDS_not_defined', 'retained_intron', 'TEC']
        self.column_names = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        self.coding_genes = {}

    def download_file(self):
        print("Looking for GTF file")
        if os.path.exists(self.gtf_file):
                print("GENCODE GTF File Found")
        else:
            print("Downloading", self.gtf_file, "version", self.version)
            url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{self.version}/gencode.v{self.version}.annotation.gtf.gz"
            output_gz_file = self.gtf_file

            print("Downloading", url)
            max_attempt = 3
            attempt = 0
            while attempt < max_attempt:
                try:
                    response = requests.get(url, stream=True, timeout=15)
                    response.raise_for_status()

                    with open(output_gz_file, 'wb') as f:
                        f.write(response.content)
                    print("Downloaded", output_gz_file)
                    break
                except requests.exceptions.Timeout:
                    attempt += 1
                    print(f"Attempt to download timed out, retrying {attempt}/{max_attempt}")
                except requests.exceptions.RequestException as e:
                    print("An error occured", e) 
                    break
            if attempt == max_attempt:
                print(f"Failed to download file after {attempt} attempts")
                sys.exit("Stopping the program.")
    
    def parse_attr(self, attr_str):
        attr = {}
        for attribute in attr_str.split(';'):
            if attribute.strip():  # Ignore empty attributes
                key, value = attribute.strip().split(' ', 1)
                value = value.strip('"')
                if key in attr:
                    attr[key].append(value)
                else:
                    attr[key] = [value]

        return attr
    
    def file_reader(self):
        gtf_df = pd.read_csv(self.gtf_file, sep='\t', comment='#', header=None, names=self.column_names)
        gtf_df = gtf_df[gtf_df['feature'].isin(['gene','transcript'])]
        print("GTF file read")
        print("Looking for coding genes")

        tot_gene = 0
        counted_genes = 0
        for index, row in gtf_df.iterrows():
            if row["feature"] == 'gene':
                tot_gene += 1
                row['attribute'] = self.parse_attr(row['attribute'])
                gene_type = row['attribute'].get('gene_type', 'N/A')[0]
                if gene_type == 'protein_coding':
                    counted_genes += 1
                    gene_id = row['attribute'].get('gene_id', 'N/A')[0].split('.')[0]
                    gene_name = row['attribute'].get('gene_name', 'N/A')[0]
                    chrom = row['chrom'][3:]
                    start, end = row['start'], row['end']
                    tag = row['attribute'].get('tag', 'N/A')
                    trans_type = row['attribute'].get('transcript_type', 'N/A')[0]
                    if trans_type not in self.unexcepted_trans_type:
                        if 'readthrough_gene' not in tag:
                            self.coding_genes[gene_id] = {"gene_id": gene_id, "gene_name":gene_name, "chrom": chrom, "start": start, "end": end, "trans_id": None, "transl_type": None, "transl_id": None, "CDS": None}


        print(f"There are {tot_gene} genes in GENCODE")
        print(f"There are {counted_genes} genes with protien coding tag")
        print(f"There are {len(self.coding_genes)} protien coding genes") 
        print("Adding translation ids and translation types")
    
        print("Adding translation types")

        for i, row in gtf_df.iterrows():
            translation_type = None
            if row['feature'] == 'transcript':
                row['attribute'] = self.parse_attr(row['attribute'])
                gene_id = row['attribute'].get('gene_id', 'N/A')[0].split('.')[0]
                if gene_id in self.coding_genes:
                    tag = row['attribute'].get('tag', 'N/A')
                    translation_type = None
                        
                    for i in tag:
                        if i == "MANE_Select":
                            translation_type = "Mane Select"
                        elif i == "Ensembl_canonical":
                            translation_type = "canonical"
                    if self.coding_genes[gene_id]['transl_type'] == None:
                        self.coding_genes[gene_id]['transl_type'] = translation_type
            
    def makeTable(self):
        print(f"Making Frame: {self.output_file}")
        final_frame = pd.DataFrame(self.coding_genes).T
        final_frame.to_excel(self.output_file) 
        print("Done")

    def run(self):
        self.download_file()
        self.file_reader()
        self.makeTable()


if __name__ == "__main__":
    processor = GENCODEProcessor()
    processor.run()
