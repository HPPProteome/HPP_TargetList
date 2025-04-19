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

df['PE'] = df['PE'].fillna(0).astype(int)
### adding MS Detection column to table
df['MS Detection'] = ""
print(df[['PE', 'UniProtKB ID', 'PeptideAtlas Category', 'Observed', 'Distinct', 'Uniquely Mapping', 'MS Detection']])

categories_and_names = [
    "No detected peptides", "No detection",
    "Identical", "Identical",
    "At least 2 unique", "2+ unique",
    "Fewer than 2 unique and fewer than 5 distinct", "Few detections",
    "Similar with few unique", "High similarity"]
categories = categories_and_names[::2]
mp_pe_values = [2, 3, 4, 5, 0]
counts = {pe: {category: 0 for category in categories} for pe in mp_pe_values}

identical_twins = {}

for index, row in df.iterrows():
    peptideatlas_category = str(row['PeptideAtlas Category'])
    uniprotkb_id = row['UniProtKB ID']

    if "identical to" in peptideatlas_category:
        identical_twins[uniprotkb_id] = True

        identical_to_id = peptideatlas_category.split(" ")[2]
        identical_twins[identical_to_id] = True

for index, row in df.iterrows():
    pe = row['PE']
    detected = row['Observed']
    distinct = row['Distinct']
    uniquely_mapping = row['Uniquely Mapping']
    uniprotkb_id = row['UniProtKB ID']

    if pe in mp_pe_values:
        if pd.isna(detected):
            counts[pe]["No detected peptides"] += 1
            df.loc[index, 'MS Detection'] = "No detections"
        elif uniprotkb_id in identical_twins:
            counts[pe]["Identical"] += 1
            df.loc[index, 'MS Detection'] = "Identical"
        elif uniquely_mapping >= 2:
            counts[pe]["At least 2 unique"] += 1
            df.loc[index, 'MS Detection'] = "2+ unique"
        elif uniquely_mapping < 2 and distinct <= 4:
            counts[pe]["Fewer than 2 unique and fewer than 5 distinct"] += 1
            df.loc[index, 'MS Detection'] = "Few detections"
        else:
            counts[pe]["Similar with few unique"] += 1
            df.loc[index, 'MS Detection'] = "High similarity"
    
    elif pe == 1:
        if pd.isna(detected):
            df.loc[index, 'MS Detection'] = "No detections"
        elif uniprotkb_id in identical_twins:
            df.loc[index, 'MS Detection'] = "Identical"
        elif uniquely_mapping >= 2:
            df.loc[index, 'MS Detection'] = "2+ unique"
        elif uniquely_mapping < 2 and distinct <= 4:
            df.loc[index, 'MS Detection'] = "Few detections"
        else:
            df.loc[index, 'MS Detection'] = "High similarity"

print("Identical twins (UniProtKB IDs):", identical_twins)

print("Final counts:\n", counts)

### bar graph
bar_heights_pe2 = []
bar_heights_pe3 = []
bar_heights_pe4 = []
bar_heights_pe5 = []
bar_heights_nope = []
for pe, categories in counts.items():
    if pe != 0:
        print(f"PE{pe}")
    elif pe == 0:
        print(f"No PE")

    for category, count in categories.items():
        print(f"    {category}: {count}")
        if pe == 2:
            bar_heights_pe2.append(int(count))
        elif pe == 3:
            bar_heights_pe3.append(int(count))
        elif pe == 4:
            bar_heights_pe4.append(int(count))
        elif pe == 5:
            bar_heights_pe5.append(int(count))
        elif pe == 0:
            bar_heights_nope.append(int(count))

fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
ax.spines[['right', 'top']].set_visible(False)

category_names = categories_and_names[1::2]
category_names = category_names[::-1]

bar_heights_pe2 = bar_heights_pe2[::-1]
bar_heights_pe3 = bar_heights_pe3[::-1]
bar_heights_pe4 = bar_heights_pe4[::-1]
bar_heights_pe5 = bar_heights_pe5[::-1]
bar_heights_nope = bar_heights_nope[::-1]
bar_heights_pe2 = np.array(bar_heights_pe2)
bar_heights_pe3 = np.array(bar_heights_pe3)
bar_heights_pe4 = np.array(bar_heights_pe4)
bar_heights_pe5 = np.array(bar_heights_pe5)
bar_heights_nope = np.array(bar_heights_nope)

print("Order of arrays on y-axis: bottom to top")
print("Categories:", category_names)
print("PE2 bar heights:", bar_heights_pe2)
print("PE3 bar heights:", bar_heights_pe3)
print("PE4 bar heights:", bar_heights_pe4)
print("PE5 bar heights:", bar_heights_pe5)
print("No PE bar heights:", bar_heights_nope)

ax.barh(category_names, bar_heights_pe2, align='center', color='#2B76B1', label="PE2")
ax.barh(category_names, bar_heights_pe3, align='center', left=bar_heights_pe2, color='#FDD83A', label="PE3")
ax.barh(category_names, bar_heights_pe4, align='center', left=(bar_heights_pe2+bar_heights_pe3), color='#0000F9', label="PE4")
ax.barh(category_names, bar_heights_pe5, align='center', left=(bar_heights_pe2+bar_heights_pe3+bar_heights_pe4), color='#870F0A', label="PE5")
max_bar_heights = ax.barh(category_names, bar_heights_nope, align='center', left=(bar_heights_pe2+bar_heights_pe3+bar_heights_pe4+bar_heights_pe5), color='Black', label="No PE")

ax.legend(fontsize=14)

ax.bar_label(max_bar_heights, padding=2, fontsize=14)

plt.xlabel("Number of Entries", fontsize=19)
plt.xticks(np.arange(0, 551, 50), fontsize=19)
plt.yticks(category_names, fontsize=14)

plt.show()

### new table with MS Detection column
print("Previewing Updated Supplementary Table 1 with MS Detection column:\n", df)
df.to_excel(excel_writer="Updated_Supplementary_Table_1.xlsx", index=False)
print("Opening Updated_Supplementary_Table_1.xlsx")
os.startfile("Updated_Supplementary_Table_1.xlsx")
print("Done")