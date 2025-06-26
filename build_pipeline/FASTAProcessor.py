#Libraries
from Bio import SeqIO
import pandas as pd
import os 
import requests
import shutil
import gzip
import subprocess
import sys
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re

class FASTAProcessor():
    def __init__(self, version):
        self.gencode_fasta = "gencode.pc_translations.fa.gz"
        self.uniprot_fasta = 'uniprot.fa'
        self.gene_file = "atlasLink.xlsx"
        self.version=version
        self.gene_file = pd.read_excel(self.gene_file)

        #Creates dictionary to store data from Gencode FASTA file
        self.ensp_dict = {}
        self.symbols_dict = {}
        self.ensg_dict = {}

        self.uniprot_dict = {}

        #Includes annotations for what each column mean
        self.secound_sheet = {"Gene ID": "ENSG number, from GENCODE GTF.",
            "Gene Symbol": "Gene Symbol, from GENCODE GTF.",
            "Chromosome": "Chromosome location, from GENCODE GTF.",
            "Start": "Nucleotide the gene starts on, from GENCODE GTF.",
            "End": "Nucleotide the gene ends on, from GENCODE GTF.",
            "Translation Type": (
                "Represents how a transcript was chosen for a single gene, from GENCODE GTF. "
                "Main Select - Indicates transcript is the single agreed-upon transcript for the protein-coding gene. "
                "Canonical - Indicates it was the 'best' transcript, selected by an individual annotation source."
            ),
            "CDS Length": "Length of gene’s amino acid sequence, from GENCODE FASTA.",
            "Protein ID & Transcirpt ID": "Identifiers taken from UniProtKB. If unavailable, they are taken from GENCODE FASTA.",
            "UniProtKB ID": "Protein Coding Gene UniProtKB ID, taken from UniProtKB.",
            "Reviewed": "Status of entry taken from UniProtKB that matches the gene’s assigned UniProtKB ID.",
            "UniProtKB Name": "Gene entry name, from UniProtKB tsv.",
            "UniProtKB Symbol": "UniProtKB ID, from UniProtKB tsv.",
            "Description": "Description of UniProtKB entry, from UniProtKB tsv.",
            "Protein Length": "Length of amino acid sequence, from UniProtKB tsv.",
            "PE": (
                "Level of protein existence (1: Evidence at protein level, 2: Evidence at transcript level, "
                "3: Inferred from homology, 4: Protein predicted, and 5: Protein Uncertain), from UniProtKB."
            ),
            "HPA Highest nTPM": "Provides the highest nTPM score for a gene, from Human Peptide Atlas.",
            "Tissues with nTPM Score Above 1": (
                "Counts the number of tissues with a nTPM score above 1.0. Out of 50 tissues, from NextProt."
            ),
            "Link Made Through": (
                "Indicates whether the UniProtKB entry was linked to the GENCODE entry through the ENSG number (gene_id), its gene symbol (gene_symbol), hand selected (Hand Selected), or matched through its GENCODE and UniProtKB sequence (Sequence)."
            ),
            "Canonical Isoform": "Notes the specific isoform sequence that is displayed for the UniProtKB ID, from UniProtKB dat file.",
            "RefSeq Identifier": "RefSeq Identifier (if availible), from UniProtKB dat file",
            "Difference in Lengths": (
                "Finds the difference in length between the UniProtKB and GENCODE entries’ amino acid sequence."
            ),

            "EC Number": "EC number of the gene, from UniProtKB tsv.",
            "nTMR": "Number of transmembrane regions, from UniProtKB tsv.",
            "Signal Peptide": "Length of a signal peptide if present, from UniProtKB tsv.",
            "PeptideAtlas Category": "Category of the protein, from PeptideAtlas.",
            "PA nObs": "Number of observed peptides, from PeptideAtlas.",
            "PA Distinct Pep": "Number of distinct peptides, from PeptideAtlas.",
            "PA Uniquely Mapping": "Number of uniquely mapping peptides, from PeptideAtlas.",
            "Hydrophobicity": "Hydrophobicity of amino acid sequence for given gene, calculated with BioPython’s ProteinAnalysis.",
            "pI": "pI number of amino acid sequence for given gene, calculated with BioPython’s ProteinAnalysis.",
            "Mass Spec ID": "Indicates whether or not the UniProtKB ID has mass spectrometry experiments associated with it, from UniProtKB dat file.",
            "3D Structure": "Indictaes whether or not information about the UniProtKB ID's 3D Structure is availible, from UniProtKB dat file.",
            "Disease Variant": "Indictaes whether or not information about the UniProtKB ID corresponding protein has one or more documented variations in its sequence that are linked to a specific disease or phenotype, from UniProtKB dat file.",
            "PPI": "Number of recorded protein on protein interactions for associated UniProtKB ID, from UniProtKB dat file."
            }

        self.secound_sheet_df = pd.DataFrame(list(self.secound_sheet.items()), columns=["Column", "Description"])


        self.noCDS = []
        self.no_uniprot = []
        self.PE5 = []


        #Build new FASTA file
        self.identifier_dict = {}
        self.seq_dict = {}
        #Columns that go in the final data frame
        self.columns_to_export = [
            'gene_id', 'gene_name', 'chrom', 'start', 'end', 'transl_type',
            'CDS', 'ENSP', 'ENST',
            'uniprot_id', 'entry_type', 'entry_name', 'gene_symbol', 'description',
            'protein length', 'evidence', 'Highest nTPM Score','Tissues with nTPM Score Above 1 (/50)', 
            'found_with', 'isoform', 'refSeq Number', 'Difference in lengths',
            'EC Number', 'Num Transmembrane Regions', 'Signal Peptide', 'PeptideAtlas Category', 'Observed', 'Distinct', 'Uniquely Mapping', 'Hydrophobicity', 
            'PI', 'Mass Spec ID', '3D-Structure', 'Disease Varient', 'PPI']






    def geneFileModifier(self):
        #Changes data from float64s to Strings
        self.gene_file['gencode_symbol'] = self.gene_file['gencode_symbol'].astype(str)
        self.gene_file['trans_id'] = self.gene_file['trans_id'].astype(str)
        self.gene_file['CDS'] = self.gene_file['CDS'].astype(str)

    def downloadGENCODE(self):
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
                
                if attempt == max_attempt:
                    print("File not downloaded, request timed out too many times.")
                    sys.exit("Exiting program")

    def downloadUniProt(self):
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

    def parseFiles(self):
        #Fills 3 dictionaries from FASTA file
        #1 that can be indexed by the ENSP number
        #1 that can be indexed with the gene symbil
        #1 that can be indexed with the ENSG number


        with gzip.open(self.gencode_fasta, "rt") as unzipped_fasta:
            for record in SeqIO.parse(unzipped_fasta, "fasta"):
                header_parts = record.description.split('|')
                ensp_number = header_parts[0]
                gene_symbol = header_parts[-2]
                cds_length = header_parts[-1]
                ensg_number = header_parts[2].split('.')[0]
                enst_number = header_parts[1]
                sequence = str(record.seq)
                self.ensp_dict[ensp_number] = {"CDS_length": cds_length, "gene_symbol": gene_symbol, "enst":enst_number, "sequence": sequence}
                if gene_symbol not in self.symbols_dict:
                    self.symbols_dict[gene_symbol] = {"CDS_lengths": [],"ensts": [],"ensps": [], "gene_symbols": [], "sequences":[]}
            
                self.symbols_dict[gene_symbol]["CDS_lengths"].append(cds_length)
                self.symbols_dict[gene_symbol]["ensts"].append(enst_number)
                self.symbols_dict[gene_symbol]["ensps"].append(ensp_number)
                self.symbols_dict[gene_symbol]['gene_symbols'].append(gene_symbol)
                self.symbols_dict[gene_symbol]['sequences'].append(sequence)    

                if ensg_number not in self.ensg_dict:
                    self.ensg_dict[ensg_number]  = {"CDS_lengths": [],"ensts": [],"ensps": [], "gene_symbols": [], "sequences":[]}

                self.ensg_dict[ensg_number]["CDS_lengths"].append(cds_length)
                self.ensg_dict[ensg_number]["ensts"].append(enst_number)
                self.ensg_dict[ensg_number]["ensps"].append(ensp_number)
                self.ensg_dict[ensg_number]['gene_symbols'].append(gene_symbol)
                self.ensg_dict[ensg_number]['sequences'].append(sequence) 

        #Parses Uniprot FASTA file
        for record in SeqIO.parse(self.uniprot_fasta, "fasta"):
            header_parts = record.description.split('|')
            id = header_parts[1]
            seq = str(record.seq)
            self.uniprot_dict[id] = seq

    def addInformation(self):
        #Initializes colunms for the DataFrame <- Not technically neccisary
        self.gene_file['Hydrophobicity'] = ''
        self.gene_file['PI'] = ''
        self.gene_file['Sequence']=''
        self.gene_file['Difference in lengths'] = ''


        sameSeq = 0
        collected = 0
        noUniProtID = 0
        noGencode = 0
        noMatch = 0
        count = 0
        incorrectENSP = 0
        strippedKey = 0
        missed = 0
        noSeq = 0

        print("Genes to investigate ENSP numbers further")
        #Checks to see if the UniProt and Gencode Sequence are the same and assignes it if it is
        for index, row in self.gene_file.iterrows():
            if isinstance(row['Key Words'], str):
                self.gene_file.at[index, 'Key Words'] = re.sub(r',\s*$', '.', row['Key Words'].replace('[', "").replace("]", "").replace(";", ", ").replace("\'", "").replace(", ,", ", "))

            else:
                self.gene_file.at[index, 'Key Words'] = ""  # Assign an empty string or handle differently


            count += 1
        #Gets sequence if UniProt and Gencode Match
            if self.ensp_dict.get(row['ENSP'], {}).get('sequence', "Hello") == self.uniprot_dict.get(row['uniprot_id'], "World"):
                self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']] #purposely doesn't use .get() to throw error if seq doesn't exist
                self.gene_file.at[index, 'CDS'] = self.ensp_dict[row['ENSP']]["CDS_length"]
                self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[row['ENSP']]["gene_symbol"]
                sameSeq += 1


        #Gets ENSP numbers for longest Sequence in Gencode if no UniProtID
            elif (row['uniprot_id'] is None or pd.isna(row['uniprot_id'])):
                        posSeq = self.ensg_dict.get(row['gene_id'], {}).get('sequences', [])
                        max_len = {"len":0, "index":0}
                        for i in posSeq:
                                if len(i) > max_len['len']:
                                        max_len['len'] = len(i)
                                        max_len['index'] = posSeq.index(i)
                        i = max_len['index']
                        self.gene_file.at[index, 'Sequence'] = self.ensg_dict[row['gene_id']]['sequences'][i]
                        self.gene_file.at[index, 'CDS'] = self.ensg_dict[row['gene_id']]["CDS_lengths"][i]
                        self.gene_file.at[index, 'gencode_symbol'] = self.ensg_dict[row['gene_id']]['gene_symbols'][i]
                        self.gene_file.at[index, 'ENST'] = self.ensg_dict[row['gene_id']]['ensts'][i]
                        self.gene_file.at[index, 'ENSP'] = self.ensg_dict[row['gene_id']]['ensps'][i]
                        self.gene_file.at[index, 'Difference in lengths'] = -int((self.ensg_dict[row['gene_id']]["CDS_lengths"][i]))
                        noUniProtID += 1

        #Takes the UniProt ENSP seq from gencode even if Sequences aren't perfect
            elif row['ENSP'] in self.ensp_dict:
                        self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']] #purposely doesn't use .get() to throw error if seq doesn't exist
                        self.gene_file.at[index, 'CDS'] = self.ensp_dict[row['ENSP']]["CDS_length"]
                        self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[row['ENSP']]["gene_symbol"]
                        self.gene_file.at[index, 'Difference in lengths'] = len(self.uniprot_dict[row['uniprot_id']]) - int(self.ensp_dict[row['ENSP']]["CDS_length"])
                        incorrectENSP += 1


        #Gets the closest length sequence from Gencode to match to UniProte 
            elif self.uniprot_dict.get(row['uniprot_id']) not in self.ensg_dict.get(row['gene_id'], {}).get('sequences', []) and self.uniprot_dict.get(row['uniprot_id']) != None and self.ensg_dict.get(row['gene_id'], {}).get('sequences', []):
                bestLenMatch = {'len':10**10, 'index':0}
                posSeq = self.ensg_dict.get(row['gene_id'], {}).get('sequences', False)

                if posSeq:
                    for i in posSeq:
                        lenDif = abs(len(self.uniprot_dict[row['uniprot_id']]) - len(i))
                        if lenDif < bestLenMatch['len']:
                            bestLenMatch['len'], bestLenMatch['index'] = lenDif, posSeq.index(i)
                    i = bestLenMatch['index']
                    self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                    self.gene_file.at[index, 'CDS'] = self.ensg_dict[row['gene_id']]["CDS_lengths"][i]
                    self.gene_file.at[index, 'gencode_symbol'] = self.ensg_dict[row['gene_id']]['gene_symbols'][i]
                    self.gene_file.at[index, 'ENST'] = self.ensg_dict[row['gene_id']]['ensts'][i]
                    self.gene_file.at[index, 'ENSP'] = self.ensg_dict[row['gene_id']]['ensps'][i]
                    self.gene_file.at[index, 'Difference in lengths'] = len(self.uniprot_dict[row['uniprot_id']]) - int((self.ensg_dict[row['gene_id']]["CDS_lengths"][i]))		
                    noMatch += 1

        #Tries and get a gencode ENSP number with the same sequence and same ENSG number in UniProt
            elif (row['ENSP'] is None or pd.isna(row['ENSP'])) and not (row['uniprot_id'] is None or pd.isna(row['uniprot_id'])):
                posSeq = self.ensg_dict.get(row['gene_id'], {}).get('sequences', False)
                if posSeq:
                    for i in range(0, len(posSeq)):
                        if posSeq[i] == self.uniprot_dict.get(row['uniprot_id']):
                            self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                            self.gene_file.at[index, 'CDS'] = self.ensg_dict[row['gene_id']]["CDS_lengths"][i]
                            self.gene_file.at[index, 'gencode_symbol'] = self.ensg_dict[row['gene_id']]['gene_symbols'][i]
                            self.gene_file.at[index, 'ENST'] = self.ensg_dict[row['gene_id']]['ensts'][i]
                            self.gene_file.at[index, 'ENSP'] = self.ensg_dict[row['gene_id']]['ensps'][i]
                            collected += 1
                            break
                else:
                    self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                    self.gene_file.at[index, 'Difference in lengths'] = "N/A"
                    noGencode += 1


        #UniProt only has a stripped ENSP number: ENSP00000494177 instead of ENSP00000494177.1, tries to make the match. Will take first .x number 
            elif "." not in row['ENSP']:
                version = [f".{i}" for i in range(1, 10)]
                for i in version:
                    posKey = row['ENSP'] + i
                    if posKey in self.ensp_dict:
                        self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                        self.gene_file.at[index, 'CDS'] = self.ensp_dict[posKey]["CDS_length"]
                        self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[posKey]['gene_symbol']
                        if self.uniprot_dict[row['uniprot_id']] != self.ensp_dict[posKey]['sequence']:
                            self.gene_file.at[index, 'Difference in lengths'] = len(self.uniprot_dict[row['uniprot_id']]) - int((self.ensp_dict[posKey]["CDS_length"]))
                        strippedKey += 1

        #Just takes the UniProtSeq and ignores ENSP num. For now it is unclear why these ones don't get caught above
            else:
                print(row['gene_id'])
                self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                missed += 1		

            if self.gene_file.at[index, 'Sequence'] == '':
                noSeq += 1

        #Calculates Hydrophobicity and pI

            sequence = ''.join([aa if aa in "ACDEFGHIKLMNPQRSTVWY" else 'L' for aa in self.gene_file.at[index, 'Sequence']]) #Removes U (selenocysteine) and other non standard char
            analysis = ProteinAnalysis(sequence)
            self.gene_file.at[index, 'Hydrophobicity'] = round(analysis.gravy(), 3)
            self.gene_file.at[index, 'PI'] = round(analysis.isoelectric_point(), 3)


        print("\n -----Stats----- \n")
        print("Number of UniProt and Gencode Sequence Matches:", sameSeq) 
        print("Number of genes with a ENSP seq to match UniProt:", collected)
        print("No UniProt ID:", noUniProtID)
        print("No Gencode:", noGencode)
        print("Incomplete matches:", noMatch)
        print("UniProt ENSP not in GENCODE:", incorrectENSP)
        print("Number of ENSP numbers with no version:", strippedKey)
        print("Genes with just a UniProtSeq (investigate further):", missed)
        print("Genes with no sequence:", noSeq)

        print("Total genes:", (missed + sameSeq + collected + noUniProtID + noGencode + noMatch + incorrectENSP + strippedKey))

        print("\nTrue number of genes:", count)

    
    def specialGeneLists(self):
        unreviewed = 0
        for index, row in self.gene_file.iterrows():
                if str(row['isoform']).startswith("ENSG"):
                    self.gene_file.at[index, 'isoform'] = None
                if row['CDS'] == "nan":
                    self.noCDS.append(row)
                    self.gene_file.at[index, 'sequence'] = 'MJA' #Indicates a lack of sequence
                if pd.isnull(row['uniprot_id']):
                    self.no_uniprot.append(row)
                if row['evidence'] == 5:
                    self.PE5.append(row)
                if row['entry_type'] == False:
                    unreviewed += 1

        
        print("\nNumber of GENCODE genes not in FASTA", len(self.noCDS))

        print("Number of GENCODE genes with no UniProt connection", len(self.no_uniprot))
        print("Number of unreviewed entries", unreviewed)



    def searchTags(self):
        self.gene_file['Mass Spec ID'] = ''
        self.gene_file['3D-Structure'] = ''
        self.gene_file['Disease Varient'] = ''


        for index, row in self.gene_file.iterrows():
    
            words = row['Key Words']
            if 'Not Present' not in words and not pd.isna(row['PPI']):
                if "3D-structure" in words:
                    self.gene_file.at[index, "3D-Structure"] = "Yes"
                else:
                    self.gene_file.at[index, "3D-Structure"] = "No"   
            
                if "Proteomicsidentification" in words:
                    self.gene_file.at[index, "Mass Spec ID"] = "Yes"
                else:
                    self.gene_file.at[index, "Mass Spec ID"] = "No" 

                if "Diseasevariant" in words:
                    self.gene_file.at[index, "Disease Varient"] = "Yes"
                else:
                    self.gene_file.at[index, "Disease Varient"] = "No" 
            
            if 'Not Present' not in words and not pd.isna(row['PPI']):
                 self.gene_file.at[index, "PPI"] = int(row['PPI'])
            else:
                 self.gene_file.at[index, "PPI"] = pd.NA


    def tableMaker(self):
        #Building dataframes
        self.searchTags()
        gene_file_selected = self.gene_file[self.columns_to_export]

        gene_file_selected = gene_file_selected.rename(columns={"gene_name":"Gene Symbol", "gene_id":"Gene ID", "transl_type": "Translation Type", "gene_symbol":"UniProtKB Symbol", "uniprot_id":"UniProtKB ID", "evidence":"PE", "CDS":"CDS Length", "found_with":"Link Made Through", "entry_name":"UniProtKB Name", "chrom":"Chromosome", "description":"Description", "isoform":"Canonical Isoform","entry_type":"Reviewed", "protein length":"Protein Length", "Difference in lengths":"Difference in Lengths",
                                                                "start":"Start", "end":"End", "ENSP":"Protein ID", "ENST":"Transcript ID", "Highest nTPM Score":"HPA Highest nTPM", "Link Made Through":"ID Link", "refSeq Number":"RefSeq Identifier",
                                                                "Num Transmembrane Regions":"nTMR", "Observed":"PA nObs", "Distinct":"PA Distinct Pep", "Uniquely Mapping":"PA Uniquely Mapping",
                                                                "PI":"pI", "3D-Structure":"3D Structure", "Disease Varient":"Disease Variant"}) 

        gene_file_selected['Reviewed'] = gene_file_selected['Reviewed'].astype(pd.BooleanDtype())
        gene_file_selected['Reviewed'] = gene_file_selected['Reviewed'].map({True: 'Reviewed', False: 'Unreviewed'})

        gene_file_selected['nTMR'] = pd.to_numeric(gene_file_selected['nTMR'], errors='coerce')
        gene_file_selected['PPI'] = pd.to_numeric(gene_file_selected['PPI'], errors='coerce').astype('Int64')

        gene_file_selected['nTMR'] = gene_file_selected['nTMR'].fillna(0).astype(int).replace(0, pd.NA)
        gene_file_selected['Signal Peptide'] = gene_file_selected['Signal Peptide'].replace("['None']", pd.NA)
        gene_file_selected['RefSeq Identifier'] = gene_file_selected['RefSeq Identifier'].apply(lambda x: "" if isinstance(x, str) and x.strip() == "None" else x)
        gene_file_selected['Signal Peptide'] = gene_file_selected['Signal Peptide'].apply(lambda x: "" if isinstance(x, str) and x.strip() == "None" else x)
        gene_file_selected['EC Number'] = gene_file_selected['EC Number'].apply(lambda x: "" if isinstance(x, str) and x.strip() == "nan" else x) 
        gene_file_selected['Difference in Lengths'] = gene_file_selected['Difference in Lengths'].apply(lambda x: "" if isinstance(x, str) and x.strip() == "N/A" else x)
        
        gene_file_selected = gene_file_selected.applymap(lambda x: x.strip() if isinstance(x, str)  else x)
        #Creates Exception files

        #No uniprot entries
        noUniprot_df = pd.DataFrame(self.no_uniprot, columns=self.columns_to_export)

        #Same UniProtID

        duplicates = gene_file_selected[gene_file_selected['UniProtKB ID'].notna() & gene_file_selected.duplicated('UniProtKB ID', keep=False)].sort_values(by='UniProtKB ID').reset_index(drop=True)
        print("Number of duplicate UniProtKB IDs", duplicates.shape[0])

        #Not in Fasta
        noCDS_df = pd.DataFrame(self.noCDS, columns=self.columns_to_export)


        #PE5 Protiens
        PE5_df = pd.DataFrame(self.PE5, columns=self.columns_to_export)
        print("Number of PE 5 proteins:", PE5_df.shape[0])
        print("Making Frames")

        #Creates table to hold sequences
        print("Main table: Supplemental_table_1.xlsx")
        output_file = "Supplemental_table_1.xlsx"
        with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
            gene_file_selected.to_excel(writer, sheet_name="Gene Table", index=False)
            self.secound_sheet_df.to_excel(writer, sheet_name="Column Annotations", index=False)



        #Extra tables
        noCDS_df.to_excel("No_Fasta_entry.xlsx", index=False)
        duplicates.to_excel("Identical_UniProtKB_entries.xlsx",index=False)
        PE5_df.to_excel("PE5_Protiens.xlsx", index=False)
        noUniprot_df.to_excel("No_UniprotKB_entries.xlsx", index=False)
        print("Tables made")



    def to_fasta(self, row):
        #Builds a fasta file by using the Uniprot ID as an idenfitier.
        #A dictionary is used to prevent repeat indentifiers, a ENSP or ENSG number is used in cases of repeats.
        #Sequences are assigned to ENSG numbers this way. Uniprot is used when a UniProtKB ID is used, GENCODE is used when a ENSG number is used
        
        if row['uniprot_id'] not in self.identifier_dict and pd.notna(row['uniprot_id']):
            line = f">{row['uniprot_id']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}|{row['description']}|{row['uniprot_id']}|{row['entry_name']}|{row['gene_name']}\n{row['Sequence']}\n"
            self.identifier_dict[row['uniprot_id']] = "Used"

        elif pd.isna(row['uniprot_id']):
            line = f">{row['ENSP']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}||||{row['gene_name']}\n{row['Sequence']}\n"

        elif pd.isna(row['ENSP']) and row['uniprot_id'] in self.identifier_dict:
            line = f">{row['gene_id']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}||||{row['gene_name']}\n{row['Sequence']}\n"

        else:
            line = f">{row['ENSP']} {row['gene_id']}|{row['ENSP']}|{len(row['Sequence'])}|{row['description']}|{row['uniprot_id']}|{row['entry_name']}|{row['gene_name']}\n{row['Sequence']}\n"

        line = line.replace('nan|','|')
        return line
    
    def buildFasta(self):
        gene_df = self.gene_file.sort_values(by='ENSP', na_position='first')
        print("Writing FASTA file")
        with open('coding_genes.fasta', 'w') as f:
            f.write(''.join(gene_df.apply(self.to_fasta, axis=1)))
        print("File written as coding_genes.fasta")

    def run(self):
        self.geneFileModifier()
        self.downloadGENCODE()
        self.downloadUniProt()
        self.parseFiles()
        self.addInformation()
        self.specialGeneLists()
        self.tableMaker()
        self.buildFasta()

if __name__ == "__main__":
    processor = FASTAProcessor(47)
    processor.run()
