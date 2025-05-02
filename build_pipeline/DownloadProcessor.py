import os 
import requests
import shutil
import sys
import gzip

class fileDownloader:
    def __init__(self, version):
        self.version = version

        self.gtf_file = "gencode.annotation.gtf.gz"
        self.gencode_fasta = "gencode.pc_translations.fa.gz"
        self.uniprot_fasta = 'uniprot.fa'

        self.atlas_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetCoreProteomeMapping/query_guest_20250416-123645-854.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_guest_20250416-123645-854&rs_page_size=1000000&output_mode=tsv"
        self.backup_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins/query_guest_20250416-172516-679.tsv?apply_action=VIEWRESULTSET&rs_set_name=query_guest_20250416-172516-679&rs_page_size=1000000&output_mode=tsv"


        self.atlas_file = "peptideAtlas.tsv"
        self.backup_atlas_file = "Secoundary_peptideAtlas.tsv"

        self.rna_file = "rna_tissue_consensus.tsv.zip"

        self.uniprot_tsv = "uniprot.tsv.gz"
        self.isoformFile = "Uniprot.dat"

    def download_GencodeGenes(self):
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

    def download_GencodeFASTA(self):
        #Downloads GENCODE FASTA file if not availible
        if os.path.exists(self.gencode_fasta):
                print("GENCODE FASTA File Found")
        else:
            
                print("Downloading gencode", self.gencode_fasta)
                url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{self.version}/gencode.v{self.version}.pc_translations.fa.gz"
                output_gz_file = self.gencode_fasta
                attempt = 0
                max_attempt = 3
                print("Downloading", url)

                while attempt < max_attempt:
                    
                    try:
                        response = requests.get(url, stream=True, timeout=10)
                        response.raise_for_status()
                        with open(output_gz_file, 'wb') as f:
                                f.write(response.content)
                        print("Downloaded", output_gz_file)
                        break
                    except requests.exceptions.Timeout:
                        attempt += 1
                        print(f"Request timed out, trying again {attempt}/{max_attempt}")

                    except requests.exceptions.RequestException as e:
                        print("File not downloaded, error:", e)
                        break
    def download_UniprotFASTA(self):
         #Downlaods UniProt FASTA file
        if not os.path.exists(self.uniprot_fasta):

            print("Downloading uniprot fasta file for needed sequences")
            url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29"
            attempt = 0
            max_attempt = 3
            
            while attempt < max_attempt:
                try:
                    response = requests.get(url, stream=True, timeout=10)
                    response.raise_for_status()
                    if response.status_code == 200:
                        with open(self.uniprot_fasta, 'w') as file:
                            file.write(response.text)
                        print("FASTA file saved as uniprot.fa")
                        break
                except requests.exceptions.Timeout:
                    attempt += 1
                    print(f"Request to download file timed out, trying again {attempt}/{max_attempt}")
                except requests.exceptions.RequestException as e:
                    print("An error occured:", e)
                    sys.exit("Exiting program")
        else:
            print("File found")

    def download_PeptideAtlas(self):
        # Check and download the primary Peptide Atlas file
        if not os.path.exists(self.atlas_file):
            attempt, max_attempts = 0, 3
            while attempt < max_attempts:
                try:
                    print(f"Downloading {self.atlas_file} from Peptide Atlas")
                    response = requests.get(self.atlas_url, stream=True, timeout=20)
                    response.raise_for_status()
                    with open(self.atlas_file, "wb") as file:
                        file.write(response.content)
                    print(f"File downloaded and saved as {self.atlas_file}")
                    break
                except requests.exceptions.Timeout:
                    attempt += 1
                    print(f"Request Timed Out, trying again {attempt}/{max_attempts}")
                except requests.exceptions.RequestException as e:
                    print(f"An error occurred: {e}")
                    sys.exit("Exiting program")
        else:
            print("Peptide Atlas file found")

        # Check and download the backup Peptide Atlas file
        if not os.path.exists(self.backup_atlas_file):
            attempt, max_attempts = 0, 3
            while attempt < max_attempts:
                try:
                    print(f"Downloading {self.backup_atlas_file} from Peptide Atlas")
                    response = requests.get(self.backup_url, stream=True, timeout=20)
                    response.raise_for_status()
                    with open(self.backup_atlas_file, "wb") as file:
                        file.write(response.content)
                    print(f"File downloaded and saved as {self.backup_atlas_file}")
                    break
                except requests.exceptions.Timeout:
                    attempt += 1
                    print(f"Request Timed Out, trying again {attempt}/{max_attempts}")
                except requests.exceptions.RequestException as e:
                    print(f"An error occurred: {e}")
                    sys.exit("Exiting program")
        else:
            print("Secondary Atlas file found")

    def download_ProtienAtlas(self):
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

    def download_Uniprot(self):
        #Checks for and downloads UniProtKB tsv file
        print("Looking for uniprot gene file")
        if os.path.exists(self.uniprot_tsv):
                print("TSV  File Found")

        else:
                print(f"Downloading {self.uniprot_tsv}")
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

    def run(self):
        self.download_GencodeFASTA()
        self.download_GencodeGenes()
        self.download_PeptideAtlas()
        self.download_ProtienAtlas()
        self.download_Uniprot()
        self.download_UniprotFASTA()
        self.downloadDat()


if __name__ == "__main__":
    processor = fileDownloader(47)
    processor.run()
    


        

