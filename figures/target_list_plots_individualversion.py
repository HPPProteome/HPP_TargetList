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

if True:
    fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')

    bars_0 = ax.bar(pe_values, num_tissues_0, label="0 Tissues", width=0.75, color='#4595D7')
    bars_1 = ax.bar(pe_values, num_tissues_1, bottom=num_tissues_0, label="1 Tissue", width=0.75, color='#FA832C')

    ax.spines[['right', 'top']].set_visible(False)

    plt.xlabel("PE Value", fontsize=21)
    plt.ylabel("Number of Entries", fontsize=21)
    plt.legend(title="Number of Tissues with nTPM score > 1 (/50)", title_fontsize=17, fontsize=17)
    plt.xticks(np.arange(0, 6, 1), fontsize=21)
    plt.yticks(np.arange(0, 1101, 100), fontsize=21)

    plt.text(x=0, y=150, s=str(grouped.loc["1", "0"]), fontsize=17, ha='center', va='center')
    plt.text(x=0, y=680, s=str(grouped.loc["1", "1"]), fontsize=17, ha='center', va='center')
    plt.text(x=1, y=130, s=str(grouped.loc["2", "0"]), fontsize=17, ha='center', va='center')
    plt.text(x=1, y=308, s=str(grouped.loc["2", "1"]), fontsize=17, ha='center', va='center')
    plt.text(x=2, y=122, s=str(grouped.loc["3", "0"]), fontsize=17, ha='center', va='center')
    plt.text(x=2, y=281, s=str(grouped.loc["3", "1"]), fontsize=17, ha='center', va='center')
    plt.text(x=3, y=20, s=str(grouped.loc["4", "0"]), fontsize=17, ha='center', va='center')
    plt.text(x=3, y=66, s=str(grouped.loc["4", "1"]), fontsize=17, ha='center', va='bottom')
    plt.text(x=4, y=18, s=str(grouped.loc["5", "0"]), fontsize=17, ha='center', va='center')
    plt.text(x=4, y=46, s=str(grouped.loc["5", "1"]), fontsize=17, ha='center', va='bottom')
    plt.text(x=5, y=-3, s=str(grouped.loc["-", "1"]), fontsize=17, ha='center', va='bottom')


### continuous value histogram of protein length

protein_lengths = df['Protein Length']

print("DataFrame Protein Lengths:\n", protein_lengths)
protein_lengths_for_plot = protein_lengths.dropna()
protein_lengths_for_plot = pd.to_numeric(protein_lengths_for_plot, errors='coerce')
protein_lengths_for_plot = protein_lengths_for_plot.astype(int)

if True:
    fig1, ax1 = plt.subplots(figsize=(10, 6), layout='constrained')

    min = 0
    max = 1000
    binsize = 10
    bins = int((max-min)/binsize)
    counts, x_floor, patches = plt.hist(protein_lengths_for_plot, bins, [min, max])

    print("Protein length counts:\n", counts)
    print("Protein length x floors:\n", x_floor)

    ax1.spines[['right', 'top']].set_visible(False)

    plt.xlabel("Protein Length", fontsize=21)
    plt.ylabel("Number of Entries", fontsize=21)
    plt.xticks(np.arange(0, 35001, 100), fontsize=21)
    plt.xlim(min, max)
    plt.yticks(np.arange(0, 601, 100), fontsize=21)

blank_protein_lengths = protein_lengths.isna().sum()
print("Number of entries with no value for protein length:", blank_protein_lengths)


### continuous value histogram of pI

pI = df['PI'].astype(float)
print("DataFrame pI:\n", pI)

if True:
    fig2, ax2 = plt.subplots(figsize=(10, 6), layout='constrained')

    min = 4
    max = 12
    binsize = 0.08
    bins = int((max-min)/binsize)
    counts, x_floor, patches = plt.hist(pI, bins, [min, max])

    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa counts:\n", counts)
    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa x floors:\n", x_floor)

    ax2.spines[['right', 'top']].set_visible(False)

    plt.xlabel("pI", fontsize=21)
    plt.ylabel("Number of Entries", fontsize=21)
    plt.xticks(np.arange(0, 13, 1), fontsize=21)
    plt.xlim(min, max)
    plt.yticks(np.arange(0, 451, 50), fontsize=21)


### continuous value histogram of hydrophobicity

hydrophobicity = df['Hydrophobicity'].astype(float)
print("DataFrame Hydrophobicity:\n", hydrophobicity)

if True:
    fig3, ax3 = plt.subplots(figsize=(10, 6), layout='constrained')

    min = -2
    max = 1.5
    binsize = 0.04
    bins = int((max-min)/binsize)
    counts, x_floor, patches = plt.hist(hydrophobicity, bins, [min, max])

    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa counts:\n", counts)
    print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa x floors:\n", x_floor)

    ax3.spines[['right', 'top']].set_visible(False)

    plt.xlabel("Hydrophobicity", fontsize=21)
    plt.ylabel("Number of Entries", fontsize=21)
    plt.xticks(np.arange(-2.5, 2.6, 0.5), fontsize=21)
    plt.xlim(min, max)
    plt.yticks(np.arange(0, 1001, 100), fontsize=21)


### discrete value histogram of number of transmembrane regions

df['Num Transmembrane Regions'] = df['Num Transmembrane Regions'].fillna(-1).astype(int)
transmembrane = df.groupby(['Num Transmembrane Regions']).size()
transmembrane = transmembrane.sort_index()

print("DataFrame Number of Entries per Number of Transmembrane Regions:\n", transmembrane)

if True:
    fig4, ax4 = plt.subplots(figsize=(10, 6), layout='constrained')

    transmembrane_for_plot = transmembrane[transmembrane.index != -1]

    num_entries_transmembrane = transmembrane_for_plot.values.tolist()
    transmembrane_for_plot = transmembrane_for_plot.index.tolist()

    bars4 = ax4.bar(transmembrane_for_plot, num_entries_transmembrane)

    ax4.spines[['right', 'top']].set_visible(False)

    plt.xlabel("Number of Transmembrane Regions", fontsize=19)
    plt.ylabel("Number of Entries", fontsize=19)
    plt.xticks(np.concatenate((np.arange(0, 10, 1), np.arange(10, 39, 2))), fontsize=19)
    plt.xlim(0, 24)
    plt.yticks(np.arange(0, 2401, 200), fontsize=19)

    ax4.bar_label(bars4, label_type='edge', fontsize=13)

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

if True:
    fig5, ax5 = plt.subplots(figsize=(10, 6), layout='constrained')

    signal_peptides_for_plot = signal_peptides[signal_peptides.index > 0]
    num_entries_signal = signal_peptides_for_plot.values.astype(int).tolist()
    signal_peptides_for_plot = signal_peptides_for_plot.index.tolist()

    ax5.bar(signal_peptides_for_plot, num_entries_signal, width=1)

    ax5.spines[['right', 'top']].set_visible(False)

    plt.xlabel("Length of Signal Peptides", fontsize=21)
    plt.ylabel("Number of Entries", fontsize=21)
    plt.xticks(np.arange(0, 71, 5), fontsize=21)
    plt.yticks(np.arange(0, 281, 40), fontsize=21)

blank_signal_peptides = signal_peptides[-1]
print("Number of entries with no value for signal peptides:", blank_signal_peptides)
no_signal_peptides_count = signal_peptides[0]
print("Number of entries with no signal peptides:", no_signal_peptides_count)



plt.show()