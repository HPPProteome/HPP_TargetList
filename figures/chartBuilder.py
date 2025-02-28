import argparse
from ChartComparer import ComparisonChart


def main():
    parser = argparse.ArgumentParser(description="Builds a Comparative chart between two gene files. Assumes the charts are in the format from code in build_pipeline (build.py)")
    parser.add_argument("--previousFile", required=True, type=str, help="This is the old file used in the comparison")
    parser.add_argument("--newFile", required=True, type=str, help="This is the new file used in the comparison")  
    args = parser.parse_args()

    processor = ComparisonChart(args.previousFile, args.newFile)
    processor.run()

if __name__ == "__main__":
    main() 
