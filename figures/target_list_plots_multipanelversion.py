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

print(df)

# setting up multipanel figure
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(21, 7), layout='constrained')


### stacked histogram of number of entries with 0 tissues and 1 tissue with nTPM score > 1 per PE value

df['PE'] = df['PE'].fillna("-1").astype(int)
df['PE'] = df['PE'].astype(str).replace("-1", "-")
print(df['PE'])
df['Tissues with nTPM Score Above 1 (/50)'] = df['Tissues with nTPM Score Above 1 (/50)'].fillna("-1").astype(int)
df['Tissues with nTPM Score Above 1 (/50)'] = df['Tissues with nTPM Score Above 1 (/50)'].astype(str).replace("-1", "-")
grouped = df.groupby(['PE', 'Tissues with nTPM Score Above 1 (/50)']).size().unstack(fill_value=0)
grouped = grouped.loc[["1", "2", "3", "4", "5", "-"]]
print("Grouped DataFrame PE vs. Tissues:\n", grouped)

pe_values = grouped.index.astype(str).tolist()
num_tissues_0 = grouped["0"].tolist()
num_tissues_1 = grouped["1"].tolist()

ax6 = axes[1, 2]

bars_0 = ax6.bar(pe_values, num_tissues_0, label="0 Tissues", width=0.75, color='#4595D7')
bars_1 = ax6.bar(pe_values, num_tissues_1, bottom=num_tissues_0, label="1 Tissue", width=0.75, color='#FA832C')

ax6.spines[['right', 'top']].set_visible(False)

ax6.set_xlabel("PE Value", fontsize=12)
ax6.set_ylabel("Number of Entries", fontsize=12)
ax6.legend(title="Number of Tissues with nTPM score > 1 (/50)", title_fontsize=10, fontsize=10)
ax6.set_xticks(np.arange(0, 6, 1))
ax6.set_yticks(np.arange(0, 1101, 100))
ax6.tick_params(axis='x', labelsize=12)
ax6.tick_params(axis='y', labelsize=12)

ax6.text(x=0, y=150, s=str(grouped.loc["1", "0"]), fontsize=10, ha='center', va='center')
ax6.text(x=0, y=680, s=str(grouped.loc["1", "1"]), fontsize=10, ha='center', va='center')
ax6.text(x=1, y=130, s=str(grouped.loc["2", "0"]), fontsize=10, ha='center', va='center')
ax6.text(x=1, y=308, s=str(grouped.loc["2", "1"]), fontsize=10, ha='center', va='center')
ax6.text(x=2, y=122, s=str(grouped.loc["3", "0"]), fontsize=10, ha='center', va='center')
ax6.text(x=2, y=281, s=str(grouped.loc["3", "1"]), fontsize=10, ha='center', va='center')
ax6.text(x=3, y=20, s=str(grouped.loc["4", "0"]), fontsize=10, ha='center', va='center')
ax6.text(x=3, y=66, s=str(grouped.loc["4", "1"]), fontsize=10, ha='center', va='bottom')
ax6.text(x=4, y=18, s=str(grouped.loc["5", "0"]), fontsize=10, ha='center', va='center')
ax6.text(x=4, y=46, s=str(grouped.loc["5", "1"]), fontsize=10, ha='center', va='bottom')
ax6.text(x=5, y=-3, s=str(grouped.loc["-", "1"]), fontsize=10, ha='center', va='bottom')


### continuous value histogram of protein length

protein_lengths = df['Protein Length']

print("DataFrame Protein Lengths:\n", protein_lengths)
protein_lengths_for_plot = protein_lengths.dropna()
protein_lengths_for_plot = pd.to_numeric(protein_lengths_for_plot, errors='coerce')
protein_lengths_for_plot = protein_lengths_for_plot.astype(int)

ax1 = axes[0, 0]

min = 0
max = 1000
binsize = 10
bins = int((max-min)/binsize)
counts, x_floor, patches = ax1.hist(protein_lengths_for_plot, bins, [min, max])

print("Protein length counts:\n", counts)
print("Protein length x floors:\n", x_floor)

ax1.spines[['right', 'top']].set_visible(False)

ax1.set_xlabel("Protein Length", fontsize=12)
ax1.set_ylabel("Number of Entries", fontsize=12)
ax1.set_xticks(np.arange(0, 35001, 100))
ax1.set_xlim(min, max)
ax1.set_yticks(np.arange(0, 601, 100))
ax1.tick_params(axis='x', labelsize=12)
ax1.tick_params(axis='y', labelsize=12)

blank_protein_lengths = protein_lengths.isna().sum()
print("Number of entries with no value for protein length:", blank_protein_lengths)


### continuous value histogram of pI

pI = df['PI'].astype(float)
print("DataFrame pI:\n", pI)

ax2 = axes[0, 1]

min = 4
max = 12
binsize = 0.08
bins = int((max-min)/binsize)
counts, x_floor, patches = ax2.hist(pI, bins, [min, max])

print("pI counts:\n", counts)
print("pI x floors:\n", x_floor)

ax2.spines[['right', 'top']].set_visible(False)

ax2.set_xlabel("pI", fontsize=12)
ax2.set_ylabel("Number of Entries", fontsize=12)
ax2.set_xticks(np.arange(0, 13, 1))
ax2.set_xlim(min, max)
ax2.set_yticks(np.arange(0, 451, 50))
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)


### continuous value histogram of hydrophobicity

hydrophobicity = df['Hydrophobicity'].astype(float)
print("DataFrame Hydrophobicity:\n", hydrophobicity)

ax3 = axes[0, 2]

min = -2
max = 1.5
binsize = 0.04
bins = int((max-min)/binsize)
counts, x_floor, patches = ax3.hist(hydrophobicity, bins, [min, max])

print("Hydrophobicity counts:\n", counts)
print("Hydrophobicity x floors:\n", x_floor)

ax3.spines[['right', 'top']].set_visible(False)

ax3.set_xlabel("Hydrophobicity", fontsize=12)
ax3.set_ylabel("Number of Entries", fontsize=12)
ax3.set_xticks(np.arange(-2.5, 2.6, 0.5))
ax3.set_xlim(min, max)
ax3.set_yticks(np.arange(0, 1001, 100))
ax3.tick_params(axis='x', labelsize=12)
ax3.tick_params(axis='y', labelsize=12)


### discrete value histogram of number of transmembrane regions

df['Num Transmembrane Regions'] = df['Num Transmembrane Regions'].fillna(-1).astype(int)
transmembrane = df.groupby(['Num Transmembrane Regions']).size()
transmembrane = transmembrane.sort_index()

print("DataFrame Number of Entries per Number of Transmembrane Regions:\n", transmembrane)

ax5 = axes[1, 1]

transmembrane_for_plot = transmembrane[transmembrane.index != -1]

num_entries_transmembrane = transmembrane_for_plot.values.tolist()
transmembrane_for_plot = transmembrane_for_plot.index.tolist()

bars5 = ax5.bar(transmembrane_for_plot, num_entries_transmembrane, rasterized=True)

ax5.spines[['right', 'top']].set_visible(False)

ax5.set_xlabel("Number of Transmembrane Regions", fontsize=11)
ax5.set_ylabel("Number of Entries", fontsize=11)
ax5.set_xticks(np.concatenate((np.arange(0, 10, 1), np.arange(10, 39, 2))))
ax5.set_xlim(0, 24)
ax5.set_yticks(np.arange(0, 2401, 200))
ax5.tick_params(axis='x', labelsize=11)
ax5.tick_params(axis='y', labelsize=11)

ax5.bar_label(bars5, label_type='edge', fontsize=7.75)

blank_transmembrane = transmembrane[-1]
print("Number of entries with no value for number of transmembrane regions:", blank_transmembrane)


### discrete value histogram of length of signal peptides

df['Signal Peptide'] = df['Signal Peptide'].fillna("-1").str.strip()
df['Signal Peptide'] = df['Signal Peptide'].replace("None", "0")
signal_peptides = df.groupby(['Signal Peptide']).size()
print("Signal peptides:", signal_peptides)

mapping = {}
for s in signal_peptides.index.astype(str).tolist():
    if ".." in s:
        try:
            mapping[s] = s.split("..")[1]
        except IndexError:
            mapping[s] = s
    else:
        mapping[s] = s
signal_peptides = signal_peptides.rename(index=mapping)
signal_peptides.index = pd.to_numeric(signal_peptides.index, errors='coerce')
signal_peptides = signal_peptides.sort_index()

print("DataFrame Number of Entries per Length of Signal Peptide:\n", signal_peptides)

ax4 = axes[1, 0]

signal_peptides_for_plot = signal_peptides[signal_peptides.index > 0]
num_entries_signal = signal_peptides_for_plot.values.astype(int).tolist()
signal_peptides_for_plot = signal_peptides_for_plot.index.tolist()

ax4.bar(signal_peptides_for_plot, num_entries_signal, width=1)

ax4.spines[['right', 'top']].set_visible(False)

ax4.set_xlabel("Length of Signal Peptides", fontsize=12)
ax4.set_ylabel("Number of Entries", fontsize=12)
ax4.set_xticks(np.arange(0, 71, 5))
ax4.set_yticks(np.arange(0, 281, 40))
ax4.tick_params(axis='x', labelsize=12)
ax4.tick_params(axis='y', labelsize=12)

blank_signal_peptides = signal_peptides[-1]
print("Number of entries with no value for signal peptides:", blank_signal_peptides)
no_signal_peptides_count = signal_peptides[0]
print("Number of entries with no signal peptides:", no_signal_peptides_count)


### plot labels
ax1.text(-180, 625, "a", fontsize=20, va='top', ha='left')
ax2.text(2.5, 470, "b", fontsize=20, va='top', ha='left')
ax3.text(-2.7, 1044, "c", fontsize=20, va='top', ha='left')
ax4.text(-13, 293, "d", fontsize=20, va='top', ha='left')
ax5.text(-4.9, 2524, "e", fontsize=20, va='top', ha='left')
ax6.text(-1.92, 1177, "f", fontsize=20, va='top', ha='left')


plt.show()