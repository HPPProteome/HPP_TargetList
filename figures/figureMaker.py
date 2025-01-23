#Libraries
import matplotlib.pyplot as plt 

#Function only displays percents higher than 3% on pie chart
def percent(num):
        if num > 3:
            return str(round(num, 1)) + '%'
        else:
            return ''



import matplotlib.pyplot as plt

# Bar chart of PeptideAtlas categories
peptideatlas_categories = {
    "N/A (entry not found)": 2,
    "Not Detected": 19,
    "Subsumed": 20,
    "Identical": 6,
    "Indistinguishable": 2,
    "Marginally Distinguished": 8,
    "Weak": 2,
    "Canonical": 17,
}

label = [key for key in peptideatlas_categories]
data = list(peptideatlas_categories.values())

# Pie chart of Suggested PE Values
Suggested_PE_Values = {
    "1": 18,
    "2": 15,
    "5": 3,
}
data_pe = list(Suggested_PE_Values.values())
label_pe = list(Suggested_PE_Values.keys())

# Pie chart of the 6 "Assessment Category" slots
assessment_categories = {
    "CPC": 21,
    "TxEOEP": 9,
    "MoFLPC": 29,
    "PPC": 12,
    "Both": 2,
    "LNC": 3,
}
data_ac = list(assessment_categories.values())
label_ac = list(assessment_categories.keys())

# Bar chart of Number of Publications
publications_count = {
    "0": 36,
    "1": 16,
    "2-4": 14,
    "5-9": 4,
    "10+": 6
}
data_pub = list(publications_count.values())
label_pub = list(publications_count.keys())

# Bar chart of HPA transcription bins
hpa_categories = {
    "N/A (entry not found)": 3,
    "Not Detected": 40,
    "Low Tissue Specificity": 4,
    "Tissue Enhanced": 6,
    "Group Enriched": 7,
    "Tissue Enriched": 16
}
data_hpa = list(hpa_categories.values())
label_hpa = list(hpa_categories.keys())

# Bar chart of AlphaFold pLDDT Scores
alphafold_plddt_scores = {
    "N/A (entry not found)": 2,
    "0-50 (Very Low)": 11,
    "50-70 (Low)": 32,
    "70-80 (High)": 8,
    "80-90 (High)": 15,
    "90-100 (Very High)": 8
}
data_plddt = list(alphafold_plddt_scores.values())
label_plddt = list(alphafold_plddt_scores.keys())

plt.figure(figsize=(21, 14))

# PeptideAtlas Categories Bar Chart
plt.subplot(2, 3, 1)
plt.barh(label, data, color="darkblue")
plt.xlabel("Protein Entries", fontsize=15)
plt.title("PeptideAtlas Categories", fontsize=25)
ax = plt.gca()
plt.text(-0.08, 1.0, "a", transform=ax.transAxes, fontsize=18, fontweight="bold")
plt.yticks(fontsize=10)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: int(x) if x.is_integer() else ""))
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# Suggested PE Values Pie Chart
plt.subplot(2, 3, 4)
wedges, texts, autotexts = plt.pie(data_pe, labels=label_pe, autopct='%.1f%%', colors=["mediumseagreen", "mediumslateblue", "brown"], startangle=90, pctdistance=0.75)
for text in texts:
    text.set_fontsize(15)
for autotext in autotexts:
    autotext.set_fontsize(15)
plt.title("Suggested PE Values", fontsize=25)
plt.yticks(fontsize=10)
ax = plt.gca()
plt.text(-0.08, 1.0, "d", transform=ax.transAxes, fontsize=18, fontweight="bold")

# Assessment Categories Pie Chart
plt.subplot(2, 3, 6)
wedges, texts, autotexts = plt.pie(data_ac, labels=label_ac, autopct=percent, colors=["lime", "olivedrab", "steelblue", "coral", "yellow", "red"], startangle=90, pctdistance=0.7)
for text in texts:
    text.set_fontsize(15)
for autotext in autotexts:
    autotext.set_fontsize(13)
ax = plt.gca()
plt.text(-0.08, 1.0, "f", transform=ax.transAxes, fontsize=18, fontweight="bold")
plt.title("Assessment Categories", fontsize=25)

# Publications Count Bar Chart
plt.subplot(2, 3, 5)
plt.bar(label_pub, data_pub, color="indigo")
plt.xlabel("Number of Publications", fontsize=15)
plt.ylabel("Protein Entries", fontsize=15)
plt.yticks(fontsize=10)
ax = plt.gca()
plt.text(-0.08, 1.0, "e", transform=ax.transAxes, fontsize=18, fontweight="bold")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# HPA Categories Bar Chart
plt.subplot(2, 3, 3)
plt.barh(label_hpa, data_hpa, color="maroon")
plt.xlabel("Protein Entries", fontsize=15)
plt.title("HPA Categories", fontsize=20)
plt.yticks(fontsize=10)
ax = plt.gca()
plt.text(-0.08, 1.0, "c", transform=ax.transAxes, fontsize=18, fontweight="bold")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# AlphaFold pLDDT Scores Bar Chart
plt.subplot(2, 3, 2)
plt.barh(label_plddt, data_plddt, color="darkgreen")
plt.xlabel("Protein Entries", fontsize=15)
plt.title("AlphaFold pLDDT Scores", fontsize=20)
ax = plt.gca()
plt.yticks(fontsize=10)
plt.text(-0.08, 1.0, "b", transform=ax.transAxes, fontsize=18, fontweight="bold")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()

# Saves the figure as an SVG file
plt.savefig("Figure 2.svg", format="svg")
print("Figure created as Figure 2.svg")

