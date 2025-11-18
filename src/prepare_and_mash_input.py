import math
import os
import random
import re
import sys
import yaml
import subprocess
from Bio import SeqIO


def sanitize_header(header, haplotype_id="1"):
    """
    Sanitizes a given header to match the PanSN format:
    - Extracts up to the first # if present and removes unwanted characters
    - Appends '#<haplotype_id>#<contig_or_scaffold_name>'
    - Ensures the header does not exceed 50 characters for downstream compatibility
    """
    # Process header, keep only valid characters (alphabets, numbers, underscores) and truncate
    parts = header.split("#")
    # Replace hyphens with underscores
    sample_name = parts[0].replace("-", "_")
    sample_name = parts[0].replace(".", "_")
    # Keep only alphabets, numbers, and underscores
    sample_name = re.sub(r'[^A-Za-z0-9_]', '', sample_name)[:30]

    # Construct PanSN-compliant header
    contig_name = haplotype_id
    sanitized_header = f"{sample_name}#{haplotype_id}#{contig_name}"

    return sanitized_header


def rename_and_compress_fasta(
        fasta_dir,
        output_file="/output/panSN_output.fasta",
        file_extensions=("fasta", "fa", "fna", "ffn", "faa", "fsa", "fas"),
        selected_files=None
):
    """
    Processes FASTA-like files in the directory by renaming headers to PanSN format,
    merging them, compressing, and returning the headers used for later reference selection.

    Args:
        fasta_dir (str): Directory containing input FASTA/FNA/FA files.
        output_file (str): Path for merged output FASTA file.
        file_extensions (tuple): Allowed file extensions for FASTA-like input.
        selected_files (list): Optional specific files to include.
    """
    sequence_count = 0
    headers = []

    if selected_files:
        missing_files = [file for file in selected_files if not os.path.exists(os.path.join(fasta_dir, file))]
        if missing_files:
            print(f"Error: The following selected files are missing: {', '.join(missing_files)}")
            exit(1)

    try:
        # Clear the output file
        with open(output_file, 'w') as outfile:
            pass

        # Process each allowed FASTA-like file
        for file_name in os.listdir(fasta_dir):
            if not any(file_name.endswith(f".{ext}") for ext in file_extensions):
                continue

            if selected_files and file_name not in selected_files:
                continue

            fasta_path = os.path.join(fasta_dir, file_name)
            with open(fasta_path, 'r') as infile, open(output_file, 'a') as outfile:
                original_header = infile.readline().strip()

                if not original_header.startswith(">"):
                    original_header = ">" + original_header

                sanitized_header = sanitize_header(original_header[1:])
                headers.append(sanitized_header)

                outfile.write(f">{sanitized_header}\n")
                for line in infile:
                    outfile.write(line)
                sequence_count += 1

        if sequence_count == 0:
            raise ValueError("No sequences were written to the output file.")

        # Downsample and compress
        downsampled_fasta = downsample_fasta(output_file, target_bases=3_000_000)
        compressed_output = bgzip_and_index(output_file)

        return compressed_output, downsampled_fasta, sequence_count, headers

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, 0, []


def downsample_fasta(fasta_file, target_bases=3_000_000):
    downsampled_file = "/output/downsampled_panSN_output.fasta"
    try:
        total_bases = sum(len(seq) for seq in SeqIO.parse(fasta_file, "fasta"))

        if total_bases > target_bases:
            subprocess.run(
                ["seqtk", "sample", "-s100", fasta_file, f"{target_bases/total_bases:.2f}"],
                stdout=open(downsampled_file, "w"),
                check=True
            )
            return downsampled_file
        else:
            print(f"Total bases ({total_bases}) is within target; using full file.")
            with open(downsampled_file, "w") as outfile, open(fasta_file, "r") as infile:
                outfile.write(infile.read())
            return downsampled_file

    except Exception as e:
        print(f"Error during downsampling: {e}")
        exit(-1)

def bgzip_and_index(fasta_file):
    gz_output_file = f"{fasta_file}.gz" if not fasta_file.endswith(".gz") else fasta_file
    try:
        if not fasta_file.endswith(".gz"):
            subprocess.run(["bgzip", "-f", fasta_file], check=True)
        subprocess.run(["samtools", "faidx", gz_output_file], check=True)
        return gz_output_file
    except subprocess.CalledProcessError as e:
        print(f"Error compressing or indexing the file: {e}")
        exit(-1)

def run_mash_triangle(fasta_file):
    try:
        result = subprocess.run(
            ["mash", "triangle", fasta_file],
            capture_output=True, text=True, check=True
        )
        mash_output = result.stdout
        if not mash_output:
            raise ValueError("Mash output is empty.")

        divergence_values = subprocess.run(
            "sed 1,1d | tr '\t' '\n' | LC_ALL=C sort -g -k 1nr | uniq | head -n 1",
            input=mash_output, capture_output=True, text=True, shell=True
        )

        max_divergence_str = divergence_values.stdout.strip()
        if not max_divergence_str:
            raise ValueError("No valid divergence found.")

        max_divergence = float(max_divergence_str)
        return math.ceil(100 - (max_divergence * 100)) - 2

    except subprocess.CalledProcessError as e:
        print(f"Error running Mash: {e}")
        exit(-1)

def prepare_fasta_and_analyze(input_dir, selected_files=None, output_yaml="/output/params.yaml"):
    """
    Prepares FASTA files for analysis, filters by selected files if provided.
    """
    merged_fasta, downsampled_fasta, record_count, headers = rename_and_compress_fasta(input_dir, selected_files=selected_files)
    
    if not merged_fasta:
        print("FASTA preparation failed.")
        return

    percent_identity = run_mash_triangle(merged_fasta)
    gzipped_fasta = bgzip_and_index(merged_fasta) if merged_fasta else None
    chosen_reference = random.choice(headers)
    results = {
        "path_to_input": gzipped_fasta,
        "downsampled_fasta": downsampled_fasta,
        "number_of_haplotypes": record_count,
        "percent_identity": percent_identity,
        "reference_header": chosen_reference
    }

    with open(output_yaml, "w") as yaml_file:
        yaml.dump(results, yaml_file, default_flow_style=False)

    print(f"Percentage identity found: {percent_identity}%")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python prepare_and_mash_input.py <input_directory> [selected_files]")
        sys.exit(1)

    input_directory = sys.argv[1]
    selected_files = sys.argv[2:] if len(sys.argv) > 2 else None
    prepare_fasta_and_analyze(input_directory, selected_files=selected_files)
