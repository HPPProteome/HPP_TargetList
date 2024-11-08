import pandas as pd

data = [
    "ENSG00000152061", "Q5R372",
    "ENSG00000179915", "Q9ULB1",
    "ENSG00000257923", "P39880",
    "ENSG00000225830", "Q03468",
    "ENSG00000186184", "P0DPB6",
    "ENSG00000021645", "Q9Y4C0",
    "ENSG00000255529", "P0CAP2",
    "ENSG00000109113", "Q9BZG1"]

genes = data[0::2]
proteins = data[1::2]

df = pd.DataFrame({'Gene ID': genes,'UniProt ID': proteins})
df.to_csv('manualFix.tsv', sep='\t', index=False)
print("TSV file of manualy selected entries created as manualFix.tsv")

