#Libraries
from Bio import SeqIO
import pandas as pd
import os 
import requests
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

        # Lists to keep track of spectial genes
        self.noCDS = []
        self.no_uniprot = []
        self.PE5 = []

        #Used to build new FASTA file
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


        #Includes annotations for what each column mean
        self.second_sheet = {"Gene ID": "ENSG number, from GENCODE GTF.",
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
            "ID Link": (
                "Indicates whether the UniProtKB entry was linked to the GENCODE entry through the ENSG number (Gene ID), its gene symbol (Gene Symbol), hand selected (Hand Selected), or matched through its GENCODE and UniProtKB sequence (Sequence)."
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

        self.second_sheet_df = pd.DataFrame(list(self.second_sheet.items()), columns=["Column", "Description"])






    def geneFileModifier(self):
        #Changes data from float64s to Strings (any nan are kept as such)
         for col in ['gencode_symbol', 'trans_id', 'CDS']:
            if col in self.gene_file.columns:
                self.gene_file[col] = self.gene_file[col].where(self.gene_file[col].isna(), self.gene_file[col].astype(str))


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
        with open(self.uniprot_fasta, 'r') as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                header_parts = record.description.split('|')
                id = header_parts[1]
                seq = str(record.seq)
                self.uniprot_dict[id] = seq
        print(f"Parsed {len(self.ensp_dict)} GENCODE ENSP records and {len(self.uniprot_dict)} UniProt records.")


    def addInformation(self):
        #Initializes colunms for the DataFrame <- Not technically neccisary
        self.gene_file['Hydrophobicity'] = pd.NA
        self.gene_file['PI'] = pd.NA
        self.gene_file['Sequence'] = ''
        self.gene_file['Difference in lengths'] = pd.NA

        #Set colunms to accept the correct types:
        self.gene_file['gencode_symbol'] = self.gene_file['gencode_symbol'].astype('string')

        # Initializes counters
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

        # Gets 3 id numbers from gene
            ens_p = row.get('ENSP')
            unip_id = row.get('uniprot_id')
            gene_id = row.get('gene_id')

        # Checks to see if the id numbers are present
            has_unip = not pd.isna(unip_id) and str(unip_id).strip() != ''
            has_ensp = not pd.isna(ens_p) and str(ens_p).strip() != ''

        #Gets sequence if UniProt and Gencode Match
            seq_gencode = None
            if has_ensp:
                seq_gencode = self.ensp_dict.get(str(ens_p), {}).get('sequence')
            seq_unip = None

            if has_unip:
                seq_unip = self.uniprot_dict.get(str(unip_id))

            if seq_gencode is not None and seq_unip is not None and seq_gencode == seq_unip:
                self.gene_file.at[index, 'Sequence'] = seq_unip
                self.gene_file.at[index, 'CDS'] = float(self.ensp_dict[str(ens_p)]["CDS_length"])
                self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[str(ens_p)]["gene_symbol"]
                sameSeq += 1
                seq_assigned = True

            else:
                seq_assigned = False


        #Gets ENSP numbers for longest Sequence in Gencode if no UniProtID
            if not seq_assigned and (not has_unip):
                if pd.notna(gene_id) and str(gene_id) in self.ensg_dict:
                    posSeq = self.ensg_dict[str(gene_id)]['sequences']
                    if posSeq:
                        lengths = [len(s) for s in posSeq]
                        i = int(np.argmax(lengths))
                        self.gene_file.at[index, 'Sequence'] = self.ensg_dict[str(gene_id)]['sequences'][i]
                        self.gene_file.at[index, 'CDS'] = float(self.ensg_dict[str(gene_id)]["CDS_lengths"][i])
                        self.gene_file.at[index, 'gencode_symbol'] = self.ensg_dict[str(gene_id)]['gene_symbols'][i]
                        self.gene_file.at[index, 'ENST'] = self.ensg_dict[str(gene_id)]['ensts'][i]
                        self.gene_file.at[index, 'ENSP'] = self.ensg_dict[str(gene_id)]['ensps'][i]
                        self.gene_file.at[index, 'Difference in lengths'] = -int(self.ensg_dict[str(gene_id)]["CDS_lengths"][i])
                        noUniProtID += 1
                        seq_assigned = True

        #Takes the UniProt ENSP seq from gencode even if Sequences aren't perfect
            elif not seq_assigned and has_ensp and str(ens_p) in self.ensp_dict:
                if has_unip:
                    self.gene_file.at[index, 'Sequence'] = seq_unip if seq_unip is not None else self.ensp_dict[str(ens_p)]['sequence']
                else:
                    self.gene_file.at[index, 'Sequence'] = self.ensp_dict[str(ens_p)]['sequence']
                self.gene_file.at[index, 'CDS'] = float(self.ensp_dict[str(ens_p)]["CDS_length"])
                self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[str(ens_p)]["gene_symbol"]
                if has_unip and seq_unip is not None:
                    try:
                        self.gene_file.at[index, 'Difference in lengths'] = len(seq_unip) - int(self.ensp_dict[str(ens_p)]["CDS_length"])
                    except Exception:
                        self.gene_file.at[index, 'Difference in lengths'] = pd.NA
                incorrectENSP += 1
                seq_assigned = True


        #Gets the closest length sequence from Gencode to match to UniProte 
            elif not seq_assigned and has_unip and pd.notna(gene_id) and str(gene_id) in self.ensg_dict:
                posSeq = self.ensg_dict[str(gene_id)]['sequences']
                if posSeq:
                    for i, s in enumerate(posSeq):
                        if s == seq_unip:
                            self.gene_file.at[index, 'Sequence'] = seq_unip
                            self.gene_file.at[index, 'CDS'] = float(self.ensg_dict[str(gene_id)]["CDS_lengths"][i])
                            self.gene_file.at[index, 'gencode_symbol'] = self.ensg_dict[str(gene_id)]['gene_symbols'][i]
                            self.gene_file.at[index, 'ENST'] = self.ensg_dict[str(gene_id)]['ensts'][i]
                            self.gene_file.at[index, 'ENSP'] = self.ensg_dict[str(gene_id)]['ensps'][i]
                            collected += 1
                            seq_assigned = True
                            break
                if not seq_assigned:
                    self.gene_file.at[index, 'Sequence'] = seq_unip if seq_unip is not None else ''
                    self.gene_file.at[index, 'Difference in lengths'] = "N/A"
                    noGencode += 1
                    seq_assigned = True

        #Tries and get a gencode ENSP number with the same sequence and same ENSG number in UniProt
            elif (row['ENSP'] is None or pd.isna(row['ENSP'])) and not (row['uniprot_id'] is None or pd.isna(row['uniprot_id'])):
                posSeq = self.ensg_dict.get(row['gene_id'], {}).get('sequences', False)
                if posSeq:
                    for i in range(0, len(posSeq)):
                        if posSeq[i] == self.uniprot_dict.get(row['uniprot_id']):
                            self.gene_file.at[index, 'Sequence'] = self.uniprot_dict[row['uniprot_id']]
                            self.gene_file.at[index, 'CDS'] = float(self.ensg_dict[row['gene_id']]["CDS_lengths"][i])
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
            elif not seq_assigned and has_ensp and "." not in str(ens_p) and has_unip:
                found = False
                for v in range(1, 11):
                    key = f"{ens_p}.{v}"
                    if key in self.ensp_dict:
                        self.gene_file.at[index, 'Sequence'] = seq_unip if seq_unip is not None else self.ensp_dict[key]['sequence']
                        self.gene_file.at[index, 'CDS'] = self.ensp_dict[key]['CDS_length']
                        self.gene_file.at[index, 'gencode_symbol'] = self.ensp_dict[key]['gene_symbol']
                        if seq_unip is not None and self.ensp_dict[key].get('sequence') is not None:
                            if seq_unip != self.ensp_dict[key]['sequence']:
                                try:
                                    self.gene_file.at[index, 'Difference in lengths'] = len(seq_unip) - int(self.ensp_dict[key]["CDS_length"])
                                except Exception:
                                    self.gene_file.at[index, 'Difference in lengths'] = pd.NA
                        strippedKey += 1
                        found = True
                        seq_assigned = True
                        break
                if not found:
                    # fallback to UniProt seq if present
                    if has_unip:
                        self.gene_file.at[index, 'Sequence'] = seq_unip
                        missed += 1
                        seq_assigned = True

            # if after all attempts still empty, keep no sequence
            else:
                noSeq += 1

        #Calculates Hydrophobicity and pI

            sequence = ''.join([aa if aa in "ACDEFGHIKLMNPQRSTVWY" else 'L' for aa in self.gene_file.at[index, 'Sequence']]) #Removes U (selenocysteine) and other non standard char
            if len(sequence) == 0:
                print(f"WARNING: No sequence for {row['gene_id']} at index {index}")
            else:
                analysis = ProteinAnalysis(sequence)
                self.gene_file.at[index, 'Hydrophobicity'] = round(analysis.gravy(), 3)
                self.gene_file.at[index, 'PI'] = round(analysis.isoelectric_point(), 3)

        #Print gene stats
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
        self.gene_file['entry_type'] = self.gene_file['entry_type'].astype(bool)
        for index, row in self.gene_file.iterrows():
            
            iso = row.get('isoform')
            if isinstance(iso, str) and iso.startswith("ENSG"):
                self.gene_file.at[index, 'isoform'] = pd.NA
            if pd.isna(row.get('CDS')):
                self.noCDS.append(row)
                self.gene_file.at[index, 'Sequence'] = 'MJA'
            if pd.isna(row.get('uniprot_id')):
                self.no_uniprot.append(row)
            if row.get('evidence') == 5:
                self.PE5.append(row)
            if row.get('entry_type') == False:
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
       
        # Set reviewed column mapping 
        if 'Reviewed' in gene_file_selected.columns:
            gene_file_selected['Reviewed'] = gene_file_selected['Reviewed'].astype(pd.BooleanDtype()).map({True: 'Reviewed', False: 'Unreviewed'})

        # numeric conversions and cleanup
        if 'nTMR' in gene_file_selected.columns:
            gene_file_selected['nTMR'] = pd.to_numeric(gene_file_selected['nTMR'], errors='coerce').fillna(0).replace(0, pd.NA).astype('Int64', errors='ignore')
        if 'PPI' in gene_file_selected.columns:
            gene_file_selected['PPI'] = pd.to_numeric(gene_file_selected['PPI'], errors='coerce').astype('Int64', errors='ignore')

        # Fix some textual "None"/"nan" placeholders
        for col in ['Signal Peptide', 'RefSeq Identifier', 'EC Number', 'Difference in Lengths']:
            if col in gene_file_selected.columns:
                gene_file_selected[col] = gene_file_selected[col].apply(lambda x: "" if (isinstance(x, str) and x.strip().lower() in ("none", "nan", "n/a")) else x)

        #Creates Exception files

        #No uniprot entries
        noUniprot_df = pd.DataFrame(self.no_uniprot, columns=self.columns_to_export)

        #Same UniProtID
        duplicates = gene_file_selected[gene_file_selected['UniProtKB ID'].notna() & gene_file_selected.duplicated('UniProtKB ID', keep=False)].sort_values(by='UniProtKB ID').reset_index(drop=True)
        print("Number of duplicate UniProtKB IDs", duplicates.shape[0])

        #Not in Fasta
        noCDS_df = pd.DataFrame(self.noCDS, columns=self.columns_to_export)


        #PE5 Proteins
        PE5_df = pd.DataFrame(self.PE5, columns=self.columns_to_export)
        print("Number of PE 5 proteins:", PE5_df.shape[0])
        print("Making Frames")

        #Creates table to hold sequences
        print("Writing Supplemental_table_1.xlsx")
        output_file = "Supplemental_table_1.xlsx"
        with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
            gene_file_selected.to_excel(writer, sheet_name="Gene Table", index=False)
            try:
                self.second_sheet_df.to_excel(writer, sheet_name="Column Annotations", index=False)
            except Exception:
                pd.DataFrame({"Column":["Column"], "Description":["Description"]}).to_excel(writer, sheet_name="Column Annotations", index=False)



        #Extra tables
        noCDS_df.to_excel("No_Fasta_entry.xlsx", index=False)
        duplicates.to_excel("Identical_UniProtKB_entries.xlsx",index=False)
        PE5_df.to_excel("PE5_Proteins.xlsx", index=False)
        noUniprot_df.to_excel("No_UniprotKB_entries.xlsx", index=False)
        print("Tables made")



    def to_fasta(self, row):
        #Builds a fasta file by using the Uniprot ID as an idenfitier.
        #A dictionary is used to prevent repeat indentifiers, a ENSP or ENSG number is used in cases of repeats.
        #Sequences are assigned to ENSG numbers this way. Uniprot is used when a UniProtKB ID is used, GENCODE is used when a ENSG number is used
        
        uid = row.get('uniprot_id')
        gene_id = row.get('gene_id', '')
        ensp = row.get('ENSP', '')
        seq = row.get('Sequence', '')
        entry_name = row.get('entry_name', '')
        gene_name = row.get('gene_name', '')

        seq = seq if isinstance(seq, str) else ''

        if pd.notna(uid) and uid not in self.identifier_dict:
            ident = str(uid)
            line = f">{ident} {gene_id}|{ensp}|{len(seq)}|{row.get('description','')}|{uid}|{entry_name}|{gene_name}\n{seq}\n"
            self.identifier_dict[ident] = True
        elif pd.isna(uid):
            line = f">{ensp} {gene_id}|{ensp}|{len(seq)}||||{gene_name}\n{seq}\n"
        elif pd.isna(ensp) and str(uid) in self.identifier_dict:
            line = f">{gene_id} {gene_id}|{ensp}|{len(seq)}||||{gene_name}\n{seq}\n"
        else:
            line = f">{ensp} {gene_id}|{ensp}|{len(seq)}|{row.get('description','')}|{uid}|{entry_name}|{gene_name}\n{seq}\n"

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
    processor = FASTAProcessor(49)
    processor.run()
