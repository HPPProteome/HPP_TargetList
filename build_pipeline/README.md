# HPP_TargetList

HPPP Proteome: Combining GENCODE and UniProtKB datasets
The goal of this pipeline is to combine the data from all protein coding genes in GENCODE with the information in UniProtKB

Install and Run:
```bash
git clone https://github.com/HPPProteome/HPP_TargetList/ 

cd HPP_TargetList/build_pipeline 

python RUN_ME.py
```

# All Scripts and their output:

# Protein_list_builder.py:
-	Downloads the newest version GENCODE’s basic gene annotations GTF file and unzips it into gencode.v47.annotation.gtf if file is not already present.
-	Outputs an excel file: coding_protiens.xlsx, which contains all 19,433 protein coding genes in the GENCODE GTF file.

# Link_to_uniprot.py:
-	Downloads all human genes from UniProtKB in a tsv file and unzips it into uniprot.tsv if file is not already present.
-	Takes in the coding_protiens.xlsx file from Protein_list_builder.py to make connections with.
-	Outputs an excel file: uniprot_output.xlsx, which contains an updated table of the 19,433 genes with all relevant UniProtKB entries.

  
# Clean_entries.py:
-	Takes in excel file: uniprot_output.xlsx from link_to_uniprot.py and removes any unnecessary UniProtKB entries.
-	Takes in manualFix.tsv to specifically choose some (current: 8) uniprot entries to match to GENCODE entries.
-	Outputs full_table.xlsx.
-	full_table.xlsx is a table of all 19,433 genes, where most UniProtKB entries have been removed to leave 1 GENCODE gene == 1 UniProtKB entry.

# Manual_fix_list_maker.py
- Builds manualFix.tsv if clean_data.py can't find the file.
- Can easily be updated to include other entries that need manual fixes.
- Includes the GENCODE ENSG number and the correct UniProt ID for a given protein coding gene.

# Update_PE.py
- Downloads rna_tissue_consensus.tsv from Human Protein Atlas.
- Takes in full_table.xlsx and adds a column called "Suggested PE".
- For every gene in full_table.xlsx, if the gene's ENSG number is present in rna_tissue_consensus.tsv, the nTPM number for any of the available tissues is greater then 1, and the gene's uniprotKB PE score is greater than 2, a 2 is added to the "Suggested PE" column.
- Outputs: updatedPE.xlsx

# Link_to_fasta:
-	Downloads the newest version of GENCODE’s Protein-coding transcript translation sequences FASTA file and unzips it into gencode.v46.pc_translations.fa if file is not present.
-	Takes in updatedPE.xlsx from update_PE.py and links the UniProtKB ENSP numbers with the GENCODE FASTA file to get GENCODE CDS length. Gene symbols are used when a ENSP is not present.
-	Outputs: Supplimental_table_1_v47.xlsx, which is the most complete table.
-	Also outputs: No_Uniprot_Entry.xlsx, which is a collection of GENCODE entries with no corrosponding UniProt ID. No_Fasta_entry.xlsx, which is a list of all GENCODE proteins in the GTF file but not in the FASTA file. Identical_UniProtKB_entries.xlsx, which is a list of all ENSG numbers that point to the same UniProt ID. PE5_Protiens.xlsx, which is a list of all proteins with a UniProtKB PE score of 5. 


# RUN_ME.py:
-	Runs all scripts in order: protein_list_builder.py > link_to_uniprot.py > clean_data.py > update_PE.py > link_to_fasta.
-	Each script includes code to automatically run the preceding script if an expected input file is missing.


