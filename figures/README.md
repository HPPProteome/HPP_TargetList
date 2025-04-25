# List of Runable Scripts:

ChartBuilder.py

figureMaker.py

target_list_plots_individualversion.py

target_list_plots_multipanelversion.py

underrepresented_protein_groups_plot.py

MS_detection_categorization_table_plot.py


# chartBuilder.py

- Runs ChartComparer.py to compare two tables created by the build_pipline to see what genes have been gained and lost between different GENCODE versions
### How to run 

```
python chartBuilder.py --previousFile [file name] --newFile [file name]

Example:
python chartBuilder.py --previousFile SupplemntryTableV46.xlsx --newFile SupplemntryTableV47.xlsx

```

### Outputs
- Outputs an SVG file named `ComparedImage.svg` that visualizes the change between the two years. **You may have edit the code to change what year is displayed**.
Example:

![image](https://github.com/user-attachments/assets/76372b1a-7f5b-4546-919f-91e1fcc7ff74)


# FigureMaker.py

# target_list_plots_individualversion.py

- Runs Supplemental_table_1.xlsx, the latest version of the GENCODE Target List found under the "lists" folder in the repository

- Outputs six figures separately, for the purpose of ease of editing and previewing, as follows:
  1. Stacked bar graph of the number of protein entries with 0 or 1 tissue(s) with an nTPM score greater than 1 out of 50
  2. Histogram of the distribution of protein lengths
  3. Histogram of the distribution of pI values
  4. Histogram of the distribution of hydrophobicity values
  5. Bar graph of the distribution of protein entries' number of transmembrane regions (if applicable)
  6. Bar graph of the distribution of signal peptide lengths (if applicable)

# target_list_plots_multipanelversion.py
- Runs Supplemental_table_1.xlsx, the latest version of the GENCODE Target List found under the "lists" folder in the repository

- Outputs the same six figures as in target_list_plots_individualversion.py, but as a single multipanel figure (2 rows, 3 columns) with the six subplots in this order:
  a) Histogram of the distribution of protein lengths.
  b) Histogram of the distribution of pI values.
  c) Histogram of the distribution of hydrophobicity values.
  d) Bar graph of the distribution of signal peptide lengths (if applicable).
  e) Bar graph of the distribution of protein entries' number of transmembrane regions (if applicable).
  f) Stacked bar graph of the number of protein entries with 0 or 1 tissue(s) with an nTPM score greater than 1 out of 50.

# underrepresented_protein_groups_plot.py
- Runs Supplemental_table_1.xlsx, the latest version of the GENCODE Target List found under the "lists" folder in the repository

- Outputs stacked bar graph of the distribution of protein entries across various protein categories (pulled from the "Description" column of Supplementary_table_1.xlsx) that are underrepresented in the PE1 category compared to missing proteins (approximately less than 90% PE1), grouped by PE value. Categories are ordered in ascending order of % PE1 with the exception of the last categories, 'Putative uncharacterized', 'Uncharacterized', 'Putative', and 'PE2-5 Other'.

# MS_detection_categorization_table_plot.py
- Runs Supplemental_table_1.xlsx, the latest version of the GENCODE Target List found under the "lists" folder in the repository

- Outputs a figure and a table summarizing PeptideAtlas' mass spectrometry evidence:
  1. Stacked bar graph of the distribution of protein entries across PeptideAtlas MS evidence categories, grouped by PE value.
  2. Updated_Supplementary_Table_1.xlsx, which adds "MS Detection" column to Supplemental_table_1.xlsx specifying the PeptideAtlas MS evidence category for each protein entry

