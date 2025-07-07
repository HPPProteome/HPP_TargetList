import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

class msDetectionCategorization():
    def __init__(self, verbose=0):
        self.gene_table = os.path.join(os.getcwd(), "Supplemental_table_1.xlsx")
        self.verbose = verbose
        self.df = pd.DataFrame()

        self.categories_and_names = [
            "No detected peptides", "No detection",
            "Identical", "Identical",
            "At least 2 unique", "2+ unique",
            "Fewer than 2 unique and fewer than 5 distinct", "Few detections",
            "Similar with few unique", "High similarity"]
        
        self.categories = self.categories_and_names[::2]
        self.mp_pe_values = [2, 3, 4, 5, 0]
        self.counts = {pe: {category: 0 for category in self.categories} for pe in self.mp_pe_values}

        self.identical_twins = {}
        self.bar_heights_pe2 = []
        self.bar_heights_pe3 = []
        self.bar_heights_pe4 = []
        self.bar_heights_pe5 = []
        self.bar_heights_nope = []

    def readTable(self):
        self.df = pd.read_excel(self.gene_table)
        self.df['PE'] = self.df['PE'].fillna(0).astype(int)

        ### adding MS Detection column to table
        self.df['MS Detection'] = ""




    def getInfo(self):
        for index, row in self.df.iterrows():
            peptideatlas_category = str(row['PeptideAtlas Category'])
            uniprotkb_id = row['UniProtKB ID']

            if "identical to" in peptideatlas_category:
                self.identical_twins[uniprotkb_id] = True

                identical_to_id = peptideatlas_category.split(" ")[2]
                self.identical_twins[identical_to_id] = True

        for index, row in self.df.iterrows():
            pe = row['PE']
            detected = row['PA nObs']
            distinct = row['PA Distinct Pep']
            uniquely_mapping = row['PA Uniquely Mapping']
            uniprotkb_id = row['UniProtKB ID']

            if pe in self.mp_pe_values:
                if pd.isna(detected):
                    self.counts[pe]["No detected peptides"] += 1
                    self.df.loc[index, 'MS Detection'] = "No detections"
                elif uniprotkb_id in self.identical_twins:
                    self.counts[pe]["Identical"] += 1
                    self.df.loc[index, 'MS Detection'] = "Identical"
                elif uniquely_mapping >= 2:
                    self.counts[pe]["At least 2 unique"] += 1
                    self.df.loc[index, 'MS Detection'] = "2+ unique"
                elif uniquely_mapping < 2 and distinct <= 4:
                    self.counts[pe]["Fewer than 2 unique and fewer than 5 distinct"] += 1
                    self.df.loc[index, 'MS Detection'] = "Few detections"
                else:
                    self.counts[pe]["Similar with few unique"] += 1
                    self.df.loc[index, 'MS Detection'] = "High similarity"
            
            elif pe == 1:
                if pd.isna(detected):
                    self.df.loc[index, 'MS Detection'] = "No detections"
                elif uniprotkb_id in self.identical_twins:
                    self.df.loc[index, 'MS Detection'] = "Identical"
                elif uniquely_mapping >= 2:
                    self.df.loc[index, 'MS Detection'] = "2+ unique"
                elif uniquely_mapping < 2 and distinct <= 4:
                    self.df.loc[index, 'MS Detection'] = "Few detections"
                else:
                    self.df.loc[index, 'MS Detection'] = "High similarity"


            if self.verbose >= 1:
                print("Final counts:\n", self.counts)

        
    def getBarHeights(self):
        for pe, categories in self.counts.items():
            if self.verbose >= 1:
                if pe != 0:
                    print(f"PE{pe}")
                elif pe == 0:
                    print(f"No PE")
            for category, count in categories.items():
                if self.verbose >= 1:
                    print(f"    {category}: {count}")
                if pe == 2:
                    self.bar_heights_pe2.append(int(count))
                elif pe == 3:
                    self.bar_heights_pe3.append(int(count))
                elif pe == 4:
                    self.bar_heights_pe4.append(int(count))
                elif pe == 5:
                    self.bar_heights_pe5.append(int(count))
                elif pe == 0:
                    self.bar_heights_nope.append(int(count))

    def buildGraph(self):
        self.getBarHeights()
        
        fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
        ax.spines[['right', 'top']].set_visible(False)

        category_names = self.categories_and_names[1::2]
        category_names = category_names[::-1]

        bar_heights_pe2 = self.bar_heights_pe2[::-1]
        bar_heights_pe3 = self.bar_heights_pe3[::-1]
        bar_heights_pe4 = self.bar_heights_pe4[::-1]
        bar_heights_pe5 = self.bar_heights_pe5[::-1]
        bar_heights_nope = self.bar_heights_nope[::-1]
        bar_heights_pe2 = np.array(bar_heights_pe2)
        bar_heights_pe3 = np.array(bar_heights_pe3)
        bar_heights_pe4 = np.array(bar_heights_pe4)
        bar_heights_pe5 = np.array(bar_heights_pe5)
        bar_heights_nope = np.array(bar_heights_nope)


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
        print("Saving chart as ms-detection-catigorization.svg")
        plt.savefig("ms-detection-catigorization.svg", format='svg', bbox_inches='tight')

        ### new table with MS Detection column
        print("Writing table: Supplementary_Table_1_MS_categorization.xlsx")
        self.df.to_excel(excel_writer="Supplementary_Table_1_MS_categorization.xlsx", index=False)
        print("Done")

    def run(self):
        self.readTable()
        self.getInfo()
        self.buildGraph()

if __name__ == "__main__":
    processor = msDetectionCategorization()
    processor.run()