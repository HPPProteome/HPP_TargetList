import os

list_files = ["coding_protiens.xlsx", "atlasLink.xlsx", "sequence_table.xlsx", "cleaned_table.xlsx", "uniprot_output.xlsx", "updatedPE.xlsx"]
downloaded_files = ["gencode.annotation.gtf.gz", "uniprot.tsv.gz", "rna_tissue_consensus.tsv.zip", "uniprot.fa", "gencode.pc_translations.fa.gz", "peptideAtlas.tsv", "Secoundary_peptideAtlas.tsv"]

print("Offering files to delete")
print("\nDeleting interm .xlsx files") 
for file in list_files:
	answer = input(f"Would you like to delete {file}? y/n ") + " "
	if answer[0].lower() == "y":
		check = os.remove(file)
		if check:
			print(f"{file} successfully removed")

print("\nDeleting Downloaded files")
for file in downloaded_files:
        answer = input(f"Would you like to delete {file}? y/n ") + " "
        if answer[0].lower() == "y":   
                check = os.remove(file)
                if check:
                        print(f"{file} successfully removed")

print("Done deleting files")
