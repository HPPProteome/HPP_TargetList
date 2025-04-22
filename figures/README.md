# List of Runable Scripts:

ChartBuilder.py

figureMaker.py

target_list_plots_individualversion.py

target_list_plots_multipanelversion.py

underrepresented_protein_groups_plot.py


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


# target_list_plots_multipanelversion.py

# underrepresented_protein_groups_plot.py
