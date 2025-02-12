import argparse
from Bio import SeqIO

def remove_duplicates(infile, outfile):
    # Dictionary to store unique sequences
    unique_sequences = {}

    # Read input file
    for record in SeqIO.parse(infile, "fasta"):
        sequence = str(record.seq)
        if sequence not in unique_sequences:
            unique_sequences[sequence] = record

    # Write the unique sequences to the output FASTA file
    with open(outfile, "w") as output_handle:
        SeqIO.write(unique_sequences.values(), output_handle, "fasta")

    print(f"Unique sequences have been saved to {outfile}")

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a FASTA file.")
    parser.add_argument('infile', type=str, help='Path to the input FASTA file')
    parser.add_argument('outfile', type=str, help='Path to the output FASTA file')

    args = parser.parse_args()
    remove_duplicates(args.infile, args.outfile)

if __name__ == "__main__":
    main()
