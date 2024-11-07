# HPP_TargetList

HPPP Proteome: Combining GENCODE and UniProtKB datasets
The goal of this pipeline is to combine the data from all protein coding genes in GENCODE with the information in UniProtKB

Install and Run:
```bash
git clone https://github.com/HPPProteome/HPP_TargetList/ 

cd HPP_TargetList/build_pipeline 

python RUN_ME.py
```
# Notes:
- Make sure that the manualFix.tsv file is in the same folder where you are running RUN_ME.py.

# All Scripts and their output:

# Protein_list_builder.py:
-	Downloads the newest version GENCODE’s basic gene annotations GTF file and unzips it into gencode.v46.annotation.gtf if file is not already present.
-	Outputs an excel file: coding_protiens.xlsx, which contains all 19,411 protein coding genes in the GENCODE GTF file.
Link_to_uniprot.py:
-	Downloads all human genes from UniProtKB in a tsv file and unzips it into uniprot.tsv if file is not already present.
-	Takes in the coding_protiens.xlsx file from Protein_list_builder.py to make connections with.
-	Outputs an excel file: uniprot_output.xlsx, which contains an updated table of the 19,411 genes with all relevant UniProtKB entries.

  
# Clean_data.py:
-	Takes in excel file: uniprot_output.xlsx from link_to_uniprot.py and removes any unnecessary UniProtKB entries.
-	Outputs two files: look_over.xlsx and full_table.xlsx.
-	look_over.xlsx is a file that holds all genes where a single UniProtKB entry could not be isolates. (Currently 29)
-	full_table.xlsx is a table of all 19,411 genes, including those in look_over.xlsx, where most UniProtKB entries have been removed to leave 1 GENCODE gene == 1 UniProtKB entry.


# Link_to_fasta:
-	Downloads the newest version of GENCODE’s Protein-coding transcript translation sequences FASTA file and unzips it into gencode.v46.pc_translations.fa if file is not present.
-	Takes in full_table.xlsx from clean_data.py and links the UniProtKB ENSP numbers with the GENCODE FASTA file to get GENCODE CDS length. Gene symbols are used when a ENSP is not present.
-	Outputs final.xlsx, which is the most complete table.
-	Also outputs: No_Uniprot_Entry.xlsx, which is a collection of 9 GENCODE entries with no corrosponding UniProt ID. No_Fasta_entry.xlsx, which is a list of all GENCODE proteins in the GTF file but not in the FASTA file. Identical_UniProtKB_entries.xlsx, which is a list of all ENSG numbers that point to the same UniProt ID.

  
# RUN_ME.py:
-	Runs all scripts in order: protein_list_builder.py > link_to_uniprot.py > clean_data.py >  link_to_fasta.
