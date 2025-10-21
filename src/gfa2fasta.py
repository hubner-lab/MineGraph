from Bio import SeqIO
import sys
import os

def parse_gfa(gfa_file):
    sequences = []
    with open(gfa_file, 'r') as file:
        for line in file:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                sequence_id = parts[1]
                sequence = parts[2]
                sequences.append((sequence_id, sequence))
    return sequences

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_gfa_file> <output_fasta_file>")
        sys.exit(1)

    input_gfa_file = sys.argv[1]
    output_fasta_file = sys.argv[2]

    sequences = parse_gfa(input_gfa_file)

    # Combine all sequences into one
    combined_sequence = ''.join(seq[1] for seq in sequences)

    # Extract file name without extension
    file_name = os.path.splitext(os.path.basename(input_gfa_file))[0]

    # Write the combined sequence to the output FASTA file with the extracted file name as header
    with open(output_fasta_file, 'w') as output_file:
        output_file.write(f">{file_name}#1\n")
        output_file.write('{}\n'.format(combined_sequence))
