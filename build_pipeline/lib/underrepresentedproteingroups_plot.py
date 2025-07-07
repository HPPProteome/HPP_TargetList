import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

class UnderrepresentedProteinGroupsChart():
    def __init__(self, verbose=0):
        self.gene_table = os.path.join(os.getcwd(), "Supplemental_table_1.xlsx")
        self.df = pd.DataFrame()
        self.verbose = verbose
        ### list of protein groups to look for and PE1-5 values
        self.protein_groups_and_names = [
            'Vomeronasal type-1 receptor', "Vomeronasal type 1 receptor",
            'lfactory receptor', "Olfactory receptor",
            'PRAME family member', "PRAME family",
            'NUT family member', "NUT family member",
            'uclear pore complex-interacting protein', "Nuclear pore complex interacting",
            'TATA-box', "TATA-box binding",
            'Speedy protein', "Speedy",
            'NBPF', "NBPF",
            'Long intergenic non-protein coding RNA', "lincRNA",
            'TAGE family member', "cTAGE family member",
            'olgin subfamily A', "Golgin subfamily A",
            'TP53', "TP53 associated",
            'Small integral membrane protein', "Small integral membrane",
            'Taste receptor', "Taste receptor",
            'eta-defensin', "Beta defensin",
            'pseudogene', "Pseudogene",
            'RNA polymerase II subunit', "RNA polymerase II subunit",
            'G antigen', "G antigen",
            'TBC1 domain family member', "TBC1 domain family",
            'biquitin carboxyl-terminal hydrolase', "Ubiquitin C-terminal hydrolase",
            'eratin-associated protein', "Keratin associated",
            'ripartite motif', "Tripartite motif containing",
            'permatogenesis', "Spermatogenesis associated",
            #'nkyrin', "Ankyrin",
            #'ransmembrane protein', "Transmembrane",
            #'ranscription factor', "Transcription factor",
            #'inc finger', "Zinc finger",

            'Putative uncharacterized', "Putative uncharacterized",
            'Uncharacterized', "Uncharacterized",
            'Putative', "Putative"
            ]
        self.protein_groups = self.protein_groups_and_names[::2]
        self.pe_values = ['1', '2', '3', '4', '5']


        ### initializing counts for each group for each PE
        self.counts = {pe: {protein_group: 0 for protein_group in self.protein_groups} for pe in self.pe_values}
        for pe in self.pe_values:
            self.counts[pe]['Other'] = 0
        self.pe1_other = 0

        ### Initialize bar height variables
        self.bar_heights_pe1 = []
        self.bar_heights_pe2 = []
        self.bar_heights_pe3 = []
        self.bar_heights_pe4 = []
        self.bar_heights_pe5 = []

    def readTable(self):
        self.df = pd.read_excel(self.gene_table)
        self.df = self.df[['Description', 'PE']]
        self.df['Description'] = self.df['Description'].fillna("-")
        self.df['PE'] = self.df['PE'].fillna("0")
        self.df['PE'] = self.df['PE'].astype(int)
   

    def getInformationFromDataFrame(self):
        ### iterating the df and counting
        for index, row in self.df.iterrows():
            description = row['Description']
            pe = str(row['PE'])

            if pe in self.pe_values: # looking for rows that are PE1-5
                found_flag = False
                for protein_group in self.protein_groups: # for each string (protein group) in the list of protein groups...
                    if protein_group in description: # ... looking for key phrases for the protein groups in the row descriptions
                        self.counts[pe][protein_group] += 1 # adding 1 to count for that PE for that protein group
                        found_flag = True
                        break # so it doesn't look for more matches with next strings in the list in the same description
                if found_flag is False and pe > '1':
                    self.counts[pe]['Other'] += 1
                elif found_flag is False and pe == '1':
                    self.pe1_other += 1

    
    def getBarHeights(self):
    ### turning counts per PE value into a list of integers for the bar heights
        for pe, groups in self.counts.items():
            if self.verbose >= 1:
                print(f"PE{pe}")
            for protein_group, count in groups.items():
                if self.verbose >= 1:
                    print(f"    {protein_group}: {count}")
                if pe == '1':
                    self.bar_heights_pe1.append(int(count))
                elif pe == '2':
                    self.bar_heights_pe2.append(int(count))
                elif pe == '3':
                    self.bar_heights_pe3.append(int(count))
                elif pe == '4':
                    self.bar_heights_pe4.append(int(count))
                elif pe == '5':
                    self.bar_heights_pe5.append(int(count))

    
    def buildChart(self):
        ### bar chart plot
        fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
        ax.spines[['right', 'top']].set_visible(False)

        group_names = self.protein_groups_and_names[1::2] + ["PE2-5 Other"]
        group_names = group_names[::-1]

        bar_heights_pe1 = self.bar_heights_pe1[::-1]
        bar_heights_pe2 = self.bar_heights_pe2[::-1]
        bar_heights_pe3 = self.bar_heights_pe3[::-1]
        bar_heights_pe4 = self.bar_heights_pe4[::-1]
        bar_heights_pe5 = self.bar_heights_pe5[::-1]
        bar_heights_pe1 = np.array(bar_heights_pe1)
        bar_heights_pe2 = np.array(bar_heights_pe2)
        bar_heights_pe3 = np.array(bar_heights_pe3)
        bar_heights_pe4 = np.array(bar_heights_pe4)
        bar_heights_pe5 = np.array(bar_heights_pe5)


        ax.barh(group_names, bar_heights_pe1, align='center', color='#3C9F38', label="PE1")
        ax.barh(group_names, bar_heights_pe2, align='center', left=bar_heights_pe1, color='#2B76B1', label="PE2")
        ax.barh(group_names, bar_heights_pe3, align='center', left=(bar_heights_pe1+bar_heights_pe2), color='#FDD83A', label="PE3")
        ax.barh(group_names, bar_heights_pe4, align='center', left=(bar_heights_pe1+bar_heights_pe2+bar_heights_pe3), color='#0000F9', label="PE4")
        max_bar_heights = ax.barh(group_names, bar_heights_pe5, align='center', left=(bar_heights_pe1+bar_heights_pe2+bar_heights_pe3+bar_heights_pe4), color='#870F0A', label="PE5")

        total_bar_heights = bar_heights_pe1 + bar_heights_pe2 + bar_heights_pe3 + bar_heights_pe4 + bar_heights_pe5

        pe1_percentages = 100 * (bar_heights_pe1 / total_bar_heights)
        pe1_percentages[0] = 100 * (self.pe1_other / (self.pe1_other + total_bar_heights[0]))
        pe1_percentages = np.round(pe1_percentages).astype(int)

        bar_labels = []
        for total, percentage in zip(total_bar_heights, pe1_percentages):
            bar_label = f"{total} ({percentage}% PE1)"
            bar_labels.append(bar_label)
        ax.bar_label(max_bar_heights, labels=bar_labels, padding=2, fontsize=11)

        plt.legend(fontsize=13, loc='center left', bbox_to_anchor=(.8, 0.5))

        plt.xlabel("Number of Entries", fontsize=17)
        plt.xticks(np.arange(0, 525, 50), fontsize=17)
        plt.yticks(group_names, fontsize=11)

        plt.savefig("underrepresented-protein-groups.svg", format='svg', bbox_inches='tight')
        print("Chart saved as underrepresented-protein-groups.svg")
        print("Done")

    def run(self):
        self.readTable()
        self.getInformationFromDataFrame()
        self.getBarHeights()
        self.buildChart()



if __name__ == "__main__":
    processor = UnderrepresentedProteinGroupsChart()
    processor.run()