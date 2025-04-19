import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import requests
import os
url = "https://raw.githubusercontent.com/HPPProteome/HPP_TargetList/main/lists/Supplemental_table_1.xlsx"
gene_file = "Supplemental_table_1.xlsx"
if not os.path.exists(gene_file):
    print("Downloading Supplemental Table 1 from GitHub")
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(gene_file, "wb") as file:
            file.write(response.content)
        print(f"File downloaded successfully and saved as '{gene_file}'")
    else:
        print(f"Failed to download file. HTTP status code: {response.status_code}")
else:
    print("Data file found")

df = pd.read_excel('Supplemental_table_1.xlsx')

df = df[['Description', 'PE']]
df['Description'] = df['Description'].fillna("-")
df['PE'] = df['PE'].fillna("0")
df['PE'] = df['PE'].astype(int)
print(df)

### list of protein groups to look for and PE1-5 values
protein_groups_and_names = [
    'Long intergenic non-protein coding RNA', "lincRNA",
    'Vomeronasal type-1 receptor', "Vomeronasal type 1 receptor",
    'lfactory receptor', "Olfactory receptor",
    'PRAME family member', "PRAME family",
    'NUT family member', "NUT family member",
    'uclear pore complex-interacting protein', "Nuclear pore complex interacting",
    'TATA-box', "TATA-box binding",
    'NBPF', "NBPF",
    'TAGE family member', "cTAGE family member",
    'Speedy protein', "Speedy",
    'olgin subfamily A', "Golgin subfamily A",
    'pseudogene', "Pseudogene",
    'TP53', "TP53 associated",
    'Small integral membrane protein', "Small integral membrane",
    'Taste receptor', "Taste receptor",
    'eta-defensin', "Beta defensin",
    'RNA polymerase II subunit', "RNA polymerase II subunit",
    'G antigen', "G antigen",
    'TBC1 domain family member', "TBC1 domain family",
    'biquitin carboxyl-terminal hydrolase', "Ubiquitin C-terminal hydrolase",
    'permatogenesis', "Spermatogenesis associated",
    'eratin-associated protein', "Keratin associated",
    'ripartite motif', "Tripartite motif containing",
    #'nkyrin', "Ankyrin",
    #'ransmembrane protein', "Transmembrane",
    #'ranscription factor', "Transcription factor",
    #'inc finger', "Zinc finger",

    'Putative uncharacterized', "Putative uncharacterized",
    'Uncharacterized', "Uncharacterized",
    'Putative', "Putative"
    ]
protein_groups = protein_groups_and_names[::2]
pe_values = ['1', '2', '3', '4', '5']

### initializing counts for each group for each PE
counts = {pe: {protein_group: 0 for protein_group in protein_groups} for pe in pe_values}
for pe in pe_values:
    counts[pe]['Other'] = 0
pe1_other = 0

### interating the df and counting
for index, row in df.iterrows():
    description = row['Description']
    pe = str(row['PE'])

    if pe in pe_values: # looking for rows that are PE1-5
        found_flag = False
        for protein_group in protein_groups: # for each string (protein group) in the list of protein groups...
            if protein_group in description: # ... looking for key phrases for the protein groups in the row descriptions
                counts[pe][protein_group] += 1 # adding 1 to count for that PE for that protein group
                found_flag = True
                break # so it doesn't look for more matches with next strings in the list in the same description
        if found_flag is False and pe > '1':
            counts[pe]['Other'] += 1
        elif found_flag is False and pe == '1':
            pe1_other += 1

### printing counts for each protein group per PE (for PE1-5)
### turning counts per PE value into a list of integers for the bar heights
bar_heights_pe1 = []
bar_heights_pe2 = []
bar_heights_pe3 = []
bar_heights_pe4 = []
bar_heights_pe5 = []
for pe, groups in counts.items():
    print(f"PE{pe}")
    for protein_group, count in groups.items():
        print(f"    {protein_group}: {count}")
        if pe == '1':
            bar_heights_pe1.append(int(count))
        elif pe == '2':
            bar_heights_pe2.append(int(count))
        elif pe == '3':
            bar_heights_pe3.append(int(count))
        elif pe == '4':
            bar_heights_pe4.append(int(count))
        elif pe == '5':
            bar_heights_pe5.append(int(count))

### bar chart plot
fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
ax.spines[['right', 'top']].set_visible(False)

group_names = protein_groups_and_names[1::2] + ["PE2-5 Other"]
group_names = group_names[::-1]

bar_heights_pe1 = bar_heights_pe1[::-1]
bar_heights_pe2 = bar_heights_pe2[::-1]
bar_heights_pe3 = bar_heights_pe3[::-1]
bar_heights_pe4 = bar_heights_pe4[::-1]
bar_heights_pe5 = bar_heights_pe5[::-1]
bar_heights_pe1 = np.array(bar_heights_pe1)
bar_heights_pe2 = np.array(bar_heights_pe2)
bar_heights_pe3 = np.array(bar_heights_pe3)
bar_heights_pe4 = np.array(bar_heights_pe4)
bar_heights_pe5 = np.array(bar_heights_pe5)

print("Order of arrays on y-axis: bottom to top")
print("Protein groups:", group_names)
print("PE1 bar heights:", bar_heights_pe1)
print("PE2 bar heights:", bar_heights_pe2)
print("PE3 bar heights:", bar_heights_pe3)
print("PE4 bar heights:", bar_heights_pe4)
print("PE5 bar heights:", bar_heights_pe5)

ax.barh(group_names, bar_heights_pe1, align='center', color='#3C9F38', label="PE1")
ax.barh(group_names, bar_heights_pe2, align='center', left=bar_heights_pe1, color='#2B76B1', label="PE2")
ax.barh(group_names, bar_heights_pe3, align='center', left=(bar_heights_pe1+bar_heights_pe2), color='#FDD83A', label="PE3")
ax.barh(group_names, bar_heights_pe4, align='center', left=(bar_heights_pe1+bar_heights_pe2+bar_heights_pe3), color='#0000F9', label="PE4")
max_bar_heights = ax.barh(group_names, bar_heights_pe5, align='center', left=(bar_heights_pe1+bar_heights_pe2+bar_heights_pe3+bar_heights_pe4), color='#870F0A', label="PE5")

total_bar_heights = bar_heights_pe1 + bar_heights_pe2 + bar_heights_pe3 + bar_heights_pe4 + bar_heights_pe5
print("Total bar heights:", total_bar_heights)
pe1_percentages = 100 * (bar_heights_pe1 / total_bar_heights)
pe1_percentages[0] = 100 * (pe1_other / (pe1_other + total_bar_heights[0]))
pe1_percentages = np.round(pe1_percentages).astype(int)
print("PE1 percentages:", pe1_percentages)
print("PE1 Other count:", pe1_other)
bar_labels = []
for total, percentage in zip(total_bar_heights, pe1_percentages):
    bar_label = f"{total} ({percentage}% PE1)"
    bar_labels.append(bar_label)
print("Bar labels:", bar_labels)
ax.bar_label(max_bar_heights, labels=bar_labels, padding=2, fontsize=11)

ax.legend(fontsize=13)

plt.xlabel("Number of Entries", fontsize=17)
plt.xticks(np.arange(0, 525, 50), fontsize=17)
plt.yticks(group_names, fontsize=11)

plt.show()