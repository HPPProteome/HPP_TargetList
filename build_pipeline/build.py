import argparse
from GENCODEProcessor import GENCODEProcessor
from UniProtProcessor import UniProtProcessor
from ManualFixProcessor import ManualList
from CleanDataProcessor import CleanDataProcessor
from ProteinAtlasProcessor import ProteinAtlasProcessor
from PeptideAtlasProcessor import PeptideAtlasProcessor
from FASTAProcessor import FASTAProcessor


def main():
    parser = argparse.ArgumentParser(description="Process GENCODE annotation files.")
    parser.add_argument("--version", default=47, type=int, help="GENCODE version to download and process (default: 47)")
    parser.add_argument("--datFile", default="UP000005640_9606.dat", type=str, help="Defines which UniProt Dat file to read from (default: UP000005640_9606.dat)")
    args = parser.parse_args()

    # Initialize and run each processor before moving to the next one
    gencode_processor = GENCODEProcessor(version=args.version)
    gencode_processor.run()

    uniprot_processor = UniProtProcessor(isoformFile=args.datFile)
    uniprot_processor.run()

    manual_fix_processor = ManualList()
    manual_fix_processor.run()

    clean_data_processor = CleanDataProcessor()
    clean_data_processor.run()

    protein_atlas_processor = ProteinAtlasProcessor()
    protein_atlas_processor.run()

    peptide_atlas_processor = PeptideAtlasProcessor()
    peptide_atlas_processor.run()

    fasta_processor = FASTAProcessor(version=args.version)
    fasta_processor.run()


if __name__ == "__main__":
    main()

