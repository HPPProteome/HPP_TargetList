from Bio import SeqIO
import pandas as pd

version = 46
fasta_file = f"gencode.v{version}.pc_translations.fa"
genes = pd.read_excel("exceptions.xlsx")
compleate_genes = pd.read_excel("complete_genes.xlsx")
fasta_entries = []

for record in SeqIO.parse(fasta_file, "fasta"):
    description = record.description.split('|')    
    gene_id = description[2].split('.')[0]
    try:
        protein_length = int(description[-1])
    except ValueError:
        continue
    
    trans_id = description[1].split('.')[0]  # Transcript ID
    transl_id = description[-2].split('.')[0]  # Translation ID
    
    fasta_entries.append((gene_id, protein_length, trans_id, transl_id))

for index, row in genes.iterrows():
    needed_length = row['protein length']
    gene_id = row['gene_id']
    
    if not pd.isna(needed_length):
        needed_length = int(needed_length)
        
        for fasta_entry in fasta_entries:
            fasta_gene_id, protein_length, trans_id, transl_id = fasta_entry
            
            if gene_id == fasta_gene_id and needed_length == protein_length:
                genes.at[index, 'CDS'] = protein_length
                genes.at[index, 'trans_id'] = trans_id
                genes.at[index, 'transl_id'] = transl_id


handled_exceptions_lst = []

not_fixed_lst = []


gene_fields = [
    "gene_id", "gene_name", "chrom", "start", "end",
    "trans_id", "transl_type", "transl_id", "CDS",
    "gencode_symbol", "uniprot_id", "reviewed",
    "entry_name", "gene_symbol", "description",
    "protein length", "entry_type", "evidence", "found_with"] 

for index, row in genes.iterrows():
        if not isinstance(row['protein length'], list):
                protein_length = row['protein length']
                cds_value = row['CDS']
                if protein_length == cds_value:
                        handled_exceptions_lst.append(row)
                else:   
                        not_fixed_lst.append(row)




not_fixed = pd.DataFrame(not_fixed_lst,columns=gene_fields)
handled_exceptions = pd.DataFrame(handled_exceptions_lst, columns=gene_fields)
diff_len_list = []
for index, row in not_fixed.iterrows():
	cds = row['CDS']
	protein_length = row['protein length']
	diff_len_list.append(abs(cds-protein_length))
not_fixed['difference in length'] = diff_len_list


print(handled_exceptions.shape)
print(not_fixed.shape)

full_genes = pd.concat([compleate_genes, handled_exceptions, not_fixed]).drop("entry_type", axis='columns')

full_genes = full_genes.sort_values(by='chrom',key=lambda chrom: chrom.str.extract(r'(\d+|X|Y|M)')[0].apply(lambda x: float(x) if x.isdigit() else 23 if x == 'X' else 24 if x == 'Y' else 25 if x == 'M' else None))

full_genes = full_genes.reset_index(drop=True)

print("Making Frame")
handled_exceptions.to_excel("handled_exceptions.xlsx")
not_fixed.to_excel("not_found.xlsx")
full_genes.iloc[:,1:].to_excel("master_list.xlsx")
