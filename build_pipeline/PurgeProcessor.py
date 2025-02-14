import os


class fileRemover():
    def __init__(self):

        self.output_files = ["atlasLink.xlsx", "sequence_table.xlsx", "Identical_UniProtKB_entries.xlsx", "No_Fasta_entry.xlsx", 
                        "cleaned_table.xlsx", "uniprot_output.xlsx", "updatedPE.xlsx", "False_ensg_over_name.xlsx", 
                        "manualFix.tsv", "protien_coding_genes.xlsx", "PE5_Protiens.xlsx", "Supplemental_table_1.xlsx", "coding_genes.fasta"]
        self.downloaded_files = ["gencode.annotation.gtf.gz", "uniprot.tsv.gz", "rna_tissue_consensus.tsv.zip", "UP000005640_9606.dat", "UP000005640_9606.dat.gz",
                            "uniprot.fa", "gencode.pc_translations.fa.gz", "peptideAtlas.tsv", "Secoundary_peptideAtlas.tsv"]

    def removeDownloads(self):
        for file in self.downloaded_files:
            if os.path.exists(file):
                os.remove(file)

    def removeOutputs(self):
        for file in self.output_files:
            if os.path.exists(file):
                os.remove(file)


    def removeAll(self):
        self.removeDownloads()
        self.removeOutputs()
