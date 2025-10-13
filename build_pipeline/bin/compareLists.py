import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

class ComparisonChart():
    def __init__(self, oldFile, newFile):
        self.oldFile = oldFile
        self.newFile = newFile

        self.new_genes = pd.read_excel(self.newFile)
        self.old_genes = pd.read_excel(self.oldFile)


        self.new_all_genes = set()
        self.new_PE1 = set()
        self.new_MP = set()
        self.new_PE5 = set()
        self.new_symbols = {}


        self.old_all_genes = set()
        self.old_PE1 = set()
        self.old_MP = set()
        self.old_PE5 = set()
        self.old_symbols = {}

        self.changed_ensg = {}
        self.num_changed = {}

    def readData(self):
        #New Table
        for index, row in self.new_genes.iterrows():
            self.new_all_genes.add(row['Gene ID'])
            self.new_symbols[row['Gene ID']] = {"PE":row['PE'], "Symbol":row['Gene Symbol']}

            if row['PE'] == 1:
                self.new_PE1.add(row['Gene ID'])

            elif row['PE'] == 5:
                self.new_PE5.add(row['Gene ID'])

            else:
                self.new_MP.add(row['Gene ID'])

        #Old Table
        for index, row in self.old_genes.iterrows():

            self.old_all_genes.add(row['Gene ID'])
            self.old_symbols[row['Gene ID']] = {"PE":row['PE'], "Symbol":row['Gene Symbol']}

            if row['PE'] == 1:
                self.old_PE1.add(row['Gene ID'])

            elif row['PE'] == 5:
                self.old_PE5.add(row['Gene ID'])

            else:
                self.old_MP.add(row['Gene ID'])
    
    def determinNameChanges(self):

        for oldKey in self.old_symbols:
            if oldKey not in self.new_symbols:
                for newKey in self.new_symbols:
                    if self.old_symbols.get(oldKey, 'old')['Symbol'] == self.new_symbols[newKey]['Symbol']:
                        self.changed_ensg[oldKey] = newKey #Old ensg : new ensg

        print("Changed Name:")
        for key in self.changed_ensg:
            print(f"{key} --> {self.changed_ensg[key]}")
            if self.old_symbols[key]['PE'] not in self.num_changed:
                self.num_changed[self.old_symbols[key]['PE']] = 0
            self.num_changed[self.old_symbols[key]['PE']] += 1
       
    


    def draw_rectangles(self):
        fig, ax = plt.subplots()
        ax.set_axis_off()
        
        rect1 = patches.Rectangle((0, 0), 2, 4, linewidth=0, edgecolor='black', facecolor='lightgrey')
        rect2 = patches.Rectangle((4, 0), 2, 4, linewidth=0, edgecolor='black', facecolor='lightgrey')
        
        
        ax.add_patch(rect1)
        ax.add_patch(rect2)


        ax.text(1, 4.2, "2024 Release", ha='center', va='center', fontsize=14, color='black')
        ax.text(5, 4.2, "2025 Release", ha='center', va='center', fontsize=14, color='black')
        
        #Total number of genes
        ax.text(1, 3.5, f"Total: {len(self.old_all_genes)}", ha='center', va='center', fontsize=13, color='black')
        ax.text(5, 3.5, f"Total: {len(self.new_all_genes)}", ha='center', va='center', fontsize=13, color='black')
        
        
    

        #To PE 1

        ax.annotate("",
                    xy=(4, 2.75), xytext=(2.1, 1.2),
                    arrowprops=dict(color='royalblue', arrowstyle="->", lw=2))
        ax.text(2.3, 1.55, f"{len(self.new_PE1 & self.old_MP)}", ha='center', va='center', rotation= 40,fontsize=10, color='royalblue')


        ax.annotate("",
                    xy=(4, 2.75), xytext=(2.1, 0.3),
                    arrowprops=dict(color='royalblue', arrowstyle="->", lw=2))
        ax.text(2.2, 0.65, f"{len(self.new_PE1 & self.old_PE5)}", ha='center', va='center', rotation= 40,fontsize=10, color='royalblue')
    

        #To PE2-4
        ax.annotate("",
                    xy=(4, 1.2), xytext=(2.1, 2.75),
                    arrowprops=dict(color='black', arrowstyle="->", lw=2))
        ax.text(2.6, 2.5, f"{len(self.new_MP & self.old_PE1)}", ha='center', va='center', rotation = -40, fontsize=10, color='black')

        ax.annotate("",
                    xy=(4, 1.2), xytext=(2.1, 0.3),
                    arrowprops=dict(color='black', arrowstyle="->", lw=2))
        ax.text(2.6, 0.7, f"{len(self.new_MP & self.old_PE5)}", ha='center', va='center', rotation = 30, fontsize=10, color='black')


        #To PE5
        ax.annotate("",
                    xy=(4, 0.3), xytext=(2.1, 2.75),
                    arrowprops=dict(color='red', arrowstyle="->", lw=2))
        ax.text(2.2, 2.25, f"{len(self.new_PE5 & self.old_PE1)}", ha='center', va='center', rotation = -47, fontsize=10, color='red')

        ax.annotate("",
                    xy=(4, 0.3), xytext=(2.1, 1.2),
                    arrowprops=dict(color='red', arrowstyle="->", lw=2))
        ax.text(2.2, 0.95, f"{len(self.new_PE5 & self.old_MP)}", ha='center', va='center', rotation = -10, fontsize=10, color='red')



    #Total PE 1
        ax.text(1, 2.75, f"{len(self.old_PE1)} PE1", ha='center', va='center', fontsize=13, color='royalblue')
        ax.text(5, 2.75, f"{len(self.new_PE1)} PE1", ha='center', va='center', fontsize=13, color='royalblue')
        
        ax.annotate("",
                    xy=(4, 2.75), xytext=(2.1, 2.75),
                    arrowprops=dict(color='royalblue', arrowstyle="->", lw=3, linestyle="dashed"))
        ax.text(3, 2.9, f"{len(self.new_PE1 & self.old_PE1)}", ha='center', va='center', fontsize=13, color='royalblue')


        #Total PE 2-4
        ax.text(1, 1.2, f"{len(self.old_MP)} MP", ha='center', va='center', fontsize=13, color='black')
        ax.text(5, 1.2, f"{len(self.new_MP)} MP", ha='center', va='center', fontsize=13, color='black')
        
        ax.annotate("",
                    xy=(4, 1.2), xytext=(2.1, 1.2),
                    arrowprops=dict(color='black', arrowstyle="->", lw=3, linestyle="dashed"))
        ax.text(3, 1.35, f"{len(self.new_MP & self.old_MP)}", ha='center', va='center', fontsize=13, color='black')


        #Total PE 5
        ax.text(1, 0.3, f"{len(self.old_PE5)} PE5", ha='center', va='center', fontsize=13, color='red')
        ax.text(5, 0.3, f"{len(self.new_PE5)} PE5", ha='center', va='center', fontsize=13, color='red')
        
        ax.annotate("",
                    xy=(4, 0.3), xytext=(2.1, 0.3),
                    arrowprops=dict(color='red', arrowstyle="->", lw=3, linestyle="dashed"))
        ax.text(3, 0.45, f"{len(self.new_PE5 & self.old_PE5)-self.num_changed.get(5, 0)}", ha='center', va='center', fontsize=13, color='red')


        #Lost Genes
        lostGenes = self.old_all_genes - self.new_all_genes
        gainedGenes = self.new_all_genes - self.old_all_genes
        
        # lost
        
        ax.text(-0.2, -0.45, f"removed", ha='center', va='center', fontsize=11, color='black')

        if len(lostGenes & self.old_PE1) >= 0:
            ax.annotate("",
                    xy=(-1, 2.75), xytext=(0,2.75),
                    arrowprops=dict(color='royalblue', arrowstyle="-", lw=1))
            ax.text(-0.45, 2.9, f"{len(lostGenes & self.old_PE1)-self.num_changed.get(1, 0)}", ha='center', va='center', fontsize=10, color='royalblue')

        if len(lostGenes & self.old_MP) >= 0:
            ax.annotate("",
                    xy=(-1, 1.2), xytext=(0, 1.2),
                    arrowprops=dict(color='black', arrowstyle="-", lw=1))
            ax.text(-0.45, 1.35, f"{len(lostGenes & self.old_MP)- self.num_changed.get(2, 0) - self.num_changed.get(3, 0) - self.num_changed.get(4, 0)}", ha='center', va='center', fontsize=10, color='black')

        if len(lostGenes & self.old_PE5) >= 0:
            ax.annotate("",
                    xy=(-1, 0.3), xytext=(0, 0.3),
                    arrowprops=dict(color='red', arrowstyle="-", lw=1))
            ax.text(-0.45, 0.45, f"{len(lostGenes & self.old_PE5)-self.num_changed.get(5, 0)}", ha='center', va='center', fontsize=10, color='red')


        ax.annotate("",
                xy=(-0.95,-0.5), xytext=(-0.95, 2.8),
                arrowprops=dict(color='black', arrowstyle="->", lw=1))

        #Gained

        ax.text(6.9, 4.2, "added", ha='center', va='center', fontsize=11, color='black')
        
        if len(gainedGenes & self.new_PE1) >= 0:
            ax.annotate("",
                    xy=(6, 2.75), xytext=(7,2.75),
                    arrowprops=dict(color='royalblue', arrowstyle="->", lw=1))
            ax.text(6.5, 2.9, f"{len(gainedGenes & self.new_PE1)-self.num_changed.get(1, 0)}", ha='center', va='center', fontsize=10, color='royalblue')

        if len(gainedGenes & self.new_MP) >= 0:
            ax.annotate("",
                    xy=(6, 1.2), xytext=(7, 1.2),
                    arrowprops=dict(color='black', arrowstyle="->", lw=1))
            ax.text(6.5, 1.35, f"{len(gainedGenes & self.new_MP)- self.num_changed.get(2, 0) - self.num_changed.get(3, 0) - self.num_changed.get(4, 0)}", ha='center', va='center', fontsize=10, color='black')

        if len(gainedGenes & self.new_PE5) >= 0:
            ax.annotate("",
                    xy=(6, 0.3), xytext=(7, 0.3),
                    arrowprops=dict(color='red', arrowstyle="->", lw=1))
            ax.text(6.5, 0.45, f"{len(gainedGenes & self.new_PE5)-self.num_changed.get(5, 0)}", ha='center', va='center', fontsize=10, color='red')

        ax.annotate("",
                xy=(6.95, 4), xytext=(6.95,0.25),
                arrowprops=dict(color='black', arrowstyle="-", lw=1))


        #Switched ENSG
        ax.text(3.0, 3.75, f"Changed identifier", ha='center', va='center', fontsize=7, color='green')

        #PE 1
        if self.num_changed.get(1, 0) > 0:
            ax.annotate("",
                xy=(4, 2.75), xytext=(2.1, 2.75),
                arrowprops=dict(
                    color='green',
                    arrowstyle="->",
                    lw=1,
                    connectionstyle="arc3,rad=-0.5"
                ))
            
            ax.text(3.05, 3.35, f"{self.num_changed.get(1, 0)}", ha='center', va='center', fontsize=9, color='green')
        
        #PE 2-4
        if (self.num_changed.get(2, 0) + self.num_changed.get(3, 0) + self.num_changed.get(4, 0)) > 0:
            ax.annotate("",
                xy=(4, 1.2), xytext=(2.1, 1.2),
                arrowprops=dict(
                    color='green',
                    arrowstyle="->",
                    lw=1,
                    connectionstyle="arc3,rad=-0.5"
                ))
            
            ax.text(3.024, 1.78, f"{self.num_changed.get(2, 0) + self.num_changed.get(3, 0) + self.num_changed.get(4, 0)}", ha='center', va='center', fontsize=9, color='green')
        #PE 5
        if self.num_changed.get(5, 0) > 0:
            ax.annotate("",
                xy=(4, 0.3), xytext=(2.1, 0.3),
                arrowprops=dict(
                    color='green',
                    arrowstyle="->",
                    lw=1,
                    connectionstyle="arc3,rad=-0.5"
                ))
            
            ax.text(2.984, 0.998, f"{self.num_changed.get(5, 0)}", ha='center', va='center', fontsize=9, color='green')



        ax.set_xlim(-1, 7)
        ax.set_ylim(-0.55, 4.5)
        ax.set_aspect(1)
        #plt.savefig('/Users/mariea/Desktop/ComparedImage.svg', format='svg', bbox_inches='tight', pad_inches=0.1) 
        plt.show()

    def run(self):
        self.readData()
        self.determinNameChanges()
        self.draw_rectangles()

if __name__ == "__main__":
    processor = ComparisonChart("/Users/mariea/Desktop/UniProt_2025_01/UniProt_2025_01.xlsx", "/Users/mariea/Desktop/Test/Supplemental_table_1.xlsx")
    processor.run()