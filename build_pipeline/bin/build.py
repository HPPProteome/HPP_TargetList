import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from DownloadProcessor import fileDownloader
from GENCODEProcessor import GENCODEProcessor
from UniProtProcessor import UniProtProcessor
from CleanDataProcessor import CleanDataProcessor
from ProteinAtlasProcessor import ProteinAtlasProcessor
from PeptideAtlasProcessor import PeptideAtlasProcessor
from FASTAProcessor import FASTAProcessor
from PurgeProcessor import fileRemover
from underrepresentedproteingroups_plot import UnderrepresentedProteinGroupsChart
from msdetectioncategorization_tableplot import msDetectionCategorization
from targetlistplots_multipanelversion import TargetListPlots


def main():
    parser = argparse.ArgumentParser(description="Combines data from GENCODE, UniProtKB, PeptideAtlas and ProteinAtlas to build Supplementry table 1")

    parser.add_argument("--version", default=48, type=int, help="GENCODE version to download and process (default: 48)")
    parser.add_argument("--build", action="store_true", default=False, help="Runs code to build FASTA file and Supplementary Table 1")
    parser.add_argument("--purge_downloads", action="store_true", default=False, help="Removes all downloaded files used to build Supplementary Table 1/FASTA file")
    parser.add_argument("--purge_output_files", action="store_true", default=False, help="Removes all interim files created by the code used to build Supplementary Table 1/FASTA file")
    parser.add_argument("--purge_all", action="store_true", default=False, help="Removes all downloaded/output files from folder")
    parser.add_argument("--build_charts", action="store_true", default=False, help="Builds extra charts about the data from generated Supplementary Table 1")


    
    args = parser.parse_args()
    
    if args.build:
       #Initializes and runs each processor before moving to the next one
        print("\n ______DOWNLOADING ALL NEEDED FILES______ \n")
        file_processor = fileDownloader(version=args.version)
        file_processor.run()

        if not os.path.exists("protien_coding_genes.xlsx"):
            print("\n ______RUNNING GENCODE PROCESSOR______ \n")
            gencode_processor = GENCODEProcessor(version=args.version)
            gencode_processor.run()
    
        if not os.path.exists("uniprot_output.xlsx"):
            print("\n ______RUNNING UNIPROT PROCESSOR______ \n")
            uniprot_processor = UniProtProcessor()
            uniprot_processor.run()


        if not os.path.exists("cleaned_table.xlsx"):
            print("\n ______RUNNING CLEAN DATA  PROCESSOR______ \n")
            clean_data_processor = CleanDataProcessor()
            clean_data_processor.run()

        if not os.path.exists("updatedPE.xlsx"):
            print("\n ______RUNNING PROTEIN ATLAS  PROCESSOR______ \n")
            protein_atlas_processor = ProteinAtlasProcessor()
            protein_atlas_processor.run()

        if not os.path.exists("atlasLink.xlsx"):
            print("\n ______RUNNING PEPTIDE ATLAS PROCESSOR______ \n")
            peptide_atlas_processor = PeptideAtlasProcessor()
            peptide_atlas_processor.run()
        

        if not os.path.exists("Supplemental_table_1.xlsx"):
            print("\n ______RUNNING FASTA PROCESSOR______ \n")
            fasta_processor = FASTAProcessor(version=args.version)
            fasta_processor.run()
    
    if args.build_charts:
        Underrep = UnderrepresentedProteinGroupsChart()
        msDetect = msDetectionCategorization()
        targetList = TargetListPlots()

        print("\n ______BUILDING UNDERREPRESENTED-PROTEIN-GROUPS CHART______ \n")
        Underrep.run()

        print("\n ______BUILDING ms DETECTION CLASSIFICATION CHART______ \n")
        msDetect.run()

        print("\n ______BUILDING 6 PANEL HISTOGRAM FIGURE______ \n")
        targetList.run()
   
    if args.purge_downloads:
        file_remover_processor = fileRemover()
        file_remover_processor.removeDownloads()

    if args.purge_output_files:
        file_remover_processor = fileRemover()
        file_remover_processor.removeOutputs()
    
    if args.purge_all:
        file_remover_processor = fileRemover()
        file_remover_processor.removeAll()
 



if __name__ == "__main__":
    main()

