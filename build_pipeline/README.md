
#THIS IS NOT UPDATED, DO NOT LOOK HERE FOR INFO

# HPP_TargetList

HPPP Proteome: Combining GENCODE and UniProtKB datasets
The goal of this pipeline is to combine the data from all protein coding genes in GENCODE with the information in UniProtKB

Install and Run:
```bash
git clone https://github.com/HPPProteome/HPP_TargetList/ 

cd HPP_TargetList/build_pipeline 

python RUN_ME.py
```
If you have a specific file you want to run the code in:
```bash
git clone https://github.com/HPPProteome/HPP_TargetList/

mkdir -p ~/folder_name

cp HPP_TargetList/build_pipeline/RUN_ME.py ~/folder_name/

cd ~/folder_name


python RUN_ME.py
```

# All Scripts and their output:

# Protein_list_builder.py:
-	Downloads the newest version GENCODE’s basic gene annotations GTF file as gencode.v47.annotation.gtf.gz if file is not already present. Program reads directly from file.

Outputs:
- coding_protiens.xlsx, which contains all 19,433 protein coding genes in the GENCODE GTF file.

# Link_to_uniprot.py:
-	Downloads all human genes from UniProtKB into  uniprot.tsv.gz file, if file is not already present.
-	Takes in the coding_protiens.xlsx file from Protein_list_builder.py to make connections.

Outputs: 
- uniprot_output.xlsx, which contains an updated table of the 19,433 genes with all relevant UniProtKB entries.

  
# Clean_entries.py:
-	Takes in excel file: uniprot_output.xlsx from link_to_uniprot.py and removes any unnecessary UniProtKB entries.
-	Takes in manualFix.tsv to specifically choose some (current: 8) uniprot entries to match to GENCODE entries.

Outputs:
-	cleaned_table.xlsx, which is a table of all 19,433 genes, where most UniProtKB entries have been removed to leave 1 GENCODE gene == 1 UniProtKB entry.
- false_ensg_over_name.xlsx, which is a file which includes a list of all ENSG numbers where we chose an unreviwed ENSG link entry over a revewied Gene Symbol link, the reveiwed UniProt links are included with the ENSG numbers. 

# Manual_fix_list_maker.py
- Builds manualFix.tsv if clean_data.py can't find the file.
- Can easily be updated to include other entries that need manual fixes.
- Includes the GENCODE ENSG number and the correct UniProt ID for a given protein-coding gene.

Outputs:
- manualFix.tsv
  
# Update_PE.py
- Downloads rna_tissue_consensus.tsv from Human Protein Atlas.
- Takes in full_table.xlsx and adds a column called "Suggested PE".
- For every gene in full_table.xlsx, if the gene's ENSG number is present in rna_tissue_consensus.tsv, the nTPM number for any of the available tissues is greater then 1, and the gene's uniprotKB PE score is greater than 2, a 2 is added to the "Suggested PE" column.

Outputs:
- updatedPE.xlsx

# PeptideAtlas_link.py
- Downloads Peptide atlas data as peptideAtlas.tsv
- Takes in updatedPE.xlsx and using ENSG numbers, links the GENCODE protein coding genes to their Peptide Atlas entry.
- Uses the Peptide Atlas information to get the genes: PeptideAtlas Category, # of Observed Peptides, # Distinct Peptides, and # Uniquely Mappin Peptides

Outputs:
- atlasLink.xlsx

# Link_to_fasta:
-	Downloads the newest version of GENCODE’s Protein-coding transcript translation sequences FASTA file as gencode.pc_translations.fa.gz and atlasLink.xlsx from peptideAtlas_link.py.
-	 Links the UniProtKB ENSP numbers with the GENCODE FASTA file to get GENCODE CDS length. Gene symbols are used when a ENSP is not present.

Outputs:
-	No_Uniprot_Entry.xlsx, which is a collection of GENCODE entries with no corrosponding UniProt ID.
-	No_Fasta_entry.xlsx, which is a list of all GENCODE proteins in the GTF file but not in the FASTA file.
-	Identical_UniProtKB_entries.xlsx, which is a list of all ENSG numbers that point to the same UniProt ID.
-	PE5_Protiens.xlsx, which is a list of all proteins with a UniProtKB PE score of 5.
-	Sequence_table.xlsx, which is the table of all genes with their information and amino acid sequences; if an amino acid sequence is not present in the GENCODE FASTA file the sequence MJA is written instead.

# Fasta_builder.py
- Takes in sequence_table.xlsx and downloads the all human gene UniProtKB FASTA file as uniprot.fa.
- For any genes in sequence_table.xlsx with the sequence MJA, the amino acid sequence is taken from the UniProtKB FASTA file.
- Builds a FASTA file with the UniProt ID as the header.
- Format: >Uniprot ID ENSG|ENSP|sequence length|Description|UniProt ID|Entry Name|Gene Symbol
  Sequence
  If the UniProt ID is unavailable or a repeat, the ENSP numver is used as the header: >ENSP ENSG|ENSP|sequence length|Description|UniProt ID|Entry Name|Gene Symbol
  Sequence
- Uses collected sequences to calculate the PI score and Hydrophobicity of each gene and appends it to the main table taken from sequences.
  
Outputs: 
- coding_genes.fasta, format discussed above
- Supplemental_table_1.xlsx, which is the most up to date table of all genes and their collected information; excludes amino acid sequences.

# RUN_ME.py:
-	Runs all scripts in order: protein_list_builder.py > link_to_uniprot.py > manual_fix_list_maker.py > clean_entries.py > update_PE.py > peptideAtlas_link.py > link_with_fasta > fasta_builder.py.
-	Each script includes code to automatically run the preceding script if an expected input file is missing.


