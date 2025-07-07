
# HPP_TargetList/build_pipeline

HPP Proteome: Combining GENCODE and UniProtKB datasets
The goal of this pipeline is to combine the data from all protein coding genes in GENCODE with the information in UniProtKB. Along with information from PeptideAtlas and ProteinAtlas.

## Starting

Clone the repo and install dependencies:
```bash
git clone https://github.com/HPPProteome/HPP_TargetList/
cd HPP_TargetList/build_pipeline 
pip install -r requirements.txt

```

Run:
```bash
# If you have a specific file you want to run the code in:
mkdir -p ~/folder_name
cd ~/folder_name

python ~/HPP_TargetList/build_pipeline/bin/build.py --build
```

## Command Line Arguments

This script accepts the following command-line arguments:

```bash
--version            GENCODE version to download and process (default: 48) Example: `--version 27`

--build              Runs code to build FASTA file and Supplementary Table 1

--purge_downloads    Removes all downloaded files used to build Supplementary Table 1/FASTA file

--purge_output_files Removes all interim files created by the code used to build Supplementary Table 1/FASTA file

--purge_all          Removes all downloaded/output files from the folder

--build_charts       Builds extra charts about the data from the generated Supplementary Table 1
```


