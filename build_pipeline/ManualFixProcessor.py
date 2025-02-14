import pandas as pd
import os

class ManualList():
    def __init__(self):
        #Lists ENSG numbers and their corrosponding UniProtKB ID
        self.data = [
            "ENSG00000152061", "Q5R372",
            "ENSG00000179915", "Q9ULB1",
            "ENSG00000257923", "P39880",
            "ENSG00000225830", "Q03468",
            "ENSG00000186184", "P0DPB6",
            "ENSG00000021645", "Q9Y4C0",
            "ENSG00000255529", "P0CAP2",
            "ENSG00000109113", "Q9BZG1"]
        
        self.output = 'manualFix.tsv'


    def createTable(self):
        if not os.path.exists(self.output):
            genes = self.data[0::2]
            proteins = self.data[1::2]
            #Creates file to hold the hand picked UniProtKB
            df = pd.DataFrame({'Gene ID': genes,'UniProt ID': proteins})
            df.to_csv(self.output, sep='\t', index=False)
            print("TSV file of manualy selected entries created as manualFix.tsv")
        else:
            print("Manual fix file already made")

    def run(self):
        self.createTable()

if __name__ == "__main__":
    processor = ManualList()
    processor.run()
