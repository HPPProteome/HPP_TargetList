
# HPP_TargetList

HPP Proteome: Combining GENCODE and UniProtKB datasets
The goal of this pipeline is to combine the data from all protein coding genes in GENCODE with the information in UniProtKB. Along with information from PeptideAtlas and ProteinAtlas.
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

cd ~/folder_name

python build.py --build
```

# All Processors and their expected output:

# GENCODEProcessor.py:

Downloads:
-	Users given version of GENCODE’s basic gene annotations GTF file as gencode.v#.annotation.gtf.gz if file is not already present.

Prints:
- How many genes are in GENCODE.
- How many genes with the protien coding tag.
- How many actual protien coding genes there are.

Outputs:
- protien_coding_genes.xlsx: A table of all protein coding genes in the GENCODE GTF file.

# UniProtProcessor.py:
Downloads:
-	Most up to date version of all human genes from UniProtKB into uniprot.tsv.gz file, if file is not already present.
-	Most up to date UniProtKB dat file (UP000005640_9606.dat.gz) as Uniprot.dat.gz

Prints:
- Number of genes from GENCODE (This should match number printed by GENCODEProcessor.py)
- List of UniProtIDs that have no listed Cannonical Isoform in the dat file
- Number of relevant reviewed UniProtKB entries found
- Number of genes with no reviewed UniProtKB entry linked to them

Outputs: 
- uniprot_output.xlsx: A table of all genes isolated by GENCODEProcessor.py, with all relevant UniProtKB entries for each gene.

  
# CleanDataProcessor.py:

Outputs:
-	cleaned_table.xlsx: A table of all isolated GENCODE genes, where all but 1 UniProtKB entries have been removed. (1 GENCODE gene == 1 UniProtKB entry.)
- false_ensg_over_name.xlsx: A table of all ENSG numbers where we chose an unreviwed ENSG link entry over a revewied Gene Symbol link.

  
# ProteinAtlasProcessor.py

Downloads:
- rna_tissue_consensus.tsv from Human Protein Atlas.

Outputs:
- updatedPE.xlsx: If a gene's ENSG number is present in rna_tissue_consensus.tsv, the nTPM number for any of the available tissues is greater then 1, and the gene's uniprotKB PE score is greater than 2, a 2 is added to the "Suggested PE" column.

# PeptideAtlasProcessor.py
Downloads:
- Peptide Atlas data as peptideAtlas.tsv
  - Peptide Atlas Core Proteom as peptideAtlas.tsv
  - Peptide Atlas all as Secoundary_peptideAtlas.tsv
 
Prints:
- Number of genes with a different UniProtKB ID compared to our table
- Number of genes found in the secondary atlas file
- Number of genes not found in either Peptide Atlas file

Outputs:
- atlasLink.xlsx: A table including prior information with Peptide Atlas information: PeptideAtlas Category, # of Observed Peptides, # Distinct Peptides, and # Uniquely Mappin Peptides.

# FASTAProcessor.py

Downloads:
- Users given version of GENCODE’s Protein-coding transcript translation sequences FASTA file as gencode.pc_translations.fa.gz
- UniProtKB's most updated FASTA file for human genes as uniprot.fa

Prints:
- Number of genes where the UniProt and Gencode Sequence match.
- Number of genes with an ENSP sequence to match UniProt.
- Number of genes with no UniProt ID.
- Number of genes with no GENCODE entry.
- Number of genes with an incomplete sequence match (UniProtKB One taken).
- Number of genes with a UniProt ENSP not in GENCODE.
- Number of ENSP numbers with no version in GENCODE.
- Number of genes with just a UniProtKB Sequence (investigate further): 1
- Number of genes with no sequence.
- Number of GENCODE genes not in FASTA.
- Number of GENCODE genes with no UniProt connection.
- Number of genes with an unreviewed UniProtKB entry.
- Number of duplicate UniProtKB IDs.
- Number of PE 5 proteins.

Outputs:
-	No_Uniprot_Entry.xlsx: A table of GENCODE entries with no corrosponding UniProt ID.
-	No_Fasta_entry.xlsx: A table GENCODE genes in the GTF file but not in the FASTA file.
-	Identical_UniProtKB_entries.xlsx: A table of all ENSG numbers that point to the same UniProt ID.
-	PE5_Protiens.xlsx: A table of all proteins with a UniProtKB PE score of 5.
- Supplemental_table_1.xlsx: Most up to date table of all genes and their collected information; excludes amino acid sequences.
- coding_genes.fasta, format discussed below
  - Format: >Uniprot ID ENSG|ENSP|sequence length|Description|UniProt ID|Entry Name|Gene Symbol
    Sequence
    If the UniProt ID is unavailable or a repeat, the ENSP numver is used as the header: >ENSP ENSG|ENSP|sequence length|Description|UniProt ID|Entry Name|Gene Symbol
    Sequence



# build.py:
-	Runs all scripts in the above order. --build
-	Allows you to pick the GENCODE version being used. --version (default = 47)
-	Delete all built tables. --purge_output_files
-	Delete all downloaded files. --purge_downloads
-	Delete all files. --purge_all



