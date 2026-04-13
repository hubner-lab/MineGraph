import argparse
import math
import os
import random
import re
import sys
import yaml
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# IUPAC + gap characters accepted as sequence data
_IUPAC = frozenset('ACGTNRYSWKMBDHVXacgtnryswkmbdhvx-.')


def _is_sequence_line(line: str) -> bool:
    """Return True if every non-space character in *line* is an IUPAC nucleotide."""
    stripped = line.strip()
    return bool(stripped) and all(c in _IUPAC for c in stripped)


def _read_sequences_robust(fasta_path: str):
    """
    Yield raw sequence strings from a FASTA file regardless of header format.

    Handles three cases:
      1. Standard    : >header\nATCG...
      2. Naked header: NAME_WITHOUT_GT\nATCG...  (first non-seq line = record separator)
      3. Headerless  : pure sequence, no header at all

    Each contiguous block of sequence lines between separators becomes one record.
    """
    current: list = []

    with open(fasta_path) as fh:
        for raw in fh:
            line = raw.rstrip()
            if not line:
                continue
            if line.startswith('>') or not _is_sequence_line(line):
                # Header line (proper or naked) — emit accumulated sequence
                if current:
                    yield ''.join(current)
                    current = []
            else:
                current.append(line)

    if current:
        yield ''.join(current)


def _sanitize_token(token: str, max_len: int) -> str:
    """Keep only [A-Za-z0-9_] and truncate."""
    token = token.replace("-", "_").replace(".", "_")
    token = re.sub(r"[^A-Za-z0-9_]", "", token)
    return token[:max_len] if max_len and max_len > 0 else token


def pansn_header(sample_name: str, haplotype_id: str, contig_name: str) -> str:
    """
    Build a PanSN-spec header using '#' delimiter:

        sample_name#haplotype_id#contig_or_scaffold_name

    Examples:
        HG002#1#ctg1234
        SAMPLE#1#0
        SAMPLE#1#1

    Notes:
      - We sanitize tokens to avoid spaces/special chars.
      - We keep the overall header reasonably short for downstream tools.
    """
    # PanSN recommends: sample#haplotype#contig
    sample = _sanitize_token(sample_name, max_len=30) or "sample"
    hap = _sanitize_token(str(haplotype_id), max_len=10) or "1"
    contig = _sanitize_token(contig_name, max_len=50) or "0"

    header = f"{sample}#{hap}#{contig}"
    # Soft cap: keep <= 80 chars (adjust if your downstream tool is stricter)
    return header[:80]


def rename_and_compress_fasta(
    fasta_dir: str,
    output_file: str = "/output/panSN_output.fasta",
    file_extensions=("fasta", "fa", "fna", "ffn", "faa", "fsa", "fas"),
    selected_files=None,
    haplotype_id: str = "1",
):
    """
    Processes FASTA-like files in a directory:
      1) Parses ALL records from each FASTA (multi-record files supported)
      2) Renames each record header to PanSN format: sample#haplotype#contig
         - sample is derived from the file name (stem)
         - contig is an index (0..N-1) per input file
      3) Merges into one FASTA
      4) bgzip + faidx
      5) Produces a downsampled FASTA for mash (seqtk)

    Returns:
      (compressed_merged_fasta_gz, downsampled_fasta, record_count, headers)
    """
    sequence_count = 0
    headers = []

    if selected_files:
        missing_files = [
            f for f in selected_files if not os.path.exists(os.path.join(fasta_dir, f))
        ]
        if missing_files:
            print(
                f"Error: The following selected files are missing: {', '.join(missing_files)}"
            )
            sys.exit(1)

    try:
        # Clear the output file
        with open(output_file, "w"):
            pass

        for file_name in sorted(os.listdir(fasta_dir)):
            if not any(file_name.endswith(f".{ext}") for ext in file_extensions):
                continue
            if selected_files and file_name not in selected_files:
                continue

            fasta_path = os.path.join(fasta_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # file stem = sample

            # Multi-record FASTA support — robust to missing/naked headers
            with open(output_file, "a") as outfile:
                for contig_idx, seq_str in enumerate(_read_sequences_robust(fasta_path)):
                    new_id = pansn_header(
                        sample_name=sample_name,
                        haplotype_id=haplotype_id,
                        contig_name=str(contig_idx),
                    )
                    headers.append(new_id)
                    record = SeqRecord(Seq(seq_str), id=new_id, name=new_id, description="")
                    SeqIO.write(record, outfile, "fasta")
                    sequence_count += 1

        if sequence_count == 0:
            raise ValueError("No sequences were written to the output file.")

        downsampled_fasta = downsample_fasta(output_file, target_bases=3_000_000)
        compressed_output = bgzip_and_index(output_file)

        return compressed_output, downsampled_fasta, sequence_count, headers

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, 0, []


def downsample_fasta(fasta_file: str, target_bases: int = 3_000_000) -> str:
    downsampled_file = "/output/downsampled_panSN_output.fasta"
    try:
        total_bases = sum(len(seq) for seq in SeqIO.parse(fasta_file, "fasta-pearson"))

        if total_bases > target_bases:
            # Sample proportionally using seqtk
            subprocess.run(
                [
                    "seqtk",
                    "sample",
                    "-s100",
                    fasta_file,
                    f"{target_bases / total_bases:.6f}",
                ],
                stdout=open(downsampled_file, "w"),
                check=True,
            )
            return downsampled_file

        print(f"Total bases ({total_bases}) is within target; using full file.")
        with open(downsampled_file, "w") as outfile, open(fasta_file, "r") as infile:
            outfile.write(infile.read())
        return downsampled_file

    except Exception as e:
        print(f"Error during downsampling: {e}")
        sys.exit(1)


def bgzip_and_index(fasta_file: str) -> str:
    gz_output_file = f"{fasta_file}.gz" if not fasta_file.endswith(".gz") else fasta_file
    try:
        if not fasta_file.endswith(".gz"):
            subprocess.run(["bgzip", "-f", fasta_file], check=True)
        subprocess.run(["samtools", "faidx", gz_output_file], check=True)
        return gz_output_file
    except subprocess.CalledProcessError as e:
        print(f"Error compressing or indexing the file: {e}")
        sys.exit(1)


def run_mash_triangle(fasta_file: str) -> int:
    """Run `mash triangle` and convert max divergence to an approximate % identity."""
    try:
        result = subprocess.run(
            ["mash", "triangle", fasta_file],
            capture_output=True,
            text=True,
            check=True,
        )
        mash_output = result.stdout
        if not mash_output:
            raise ValueError("Mash output is empty.")

        divergence_values = subprocess.run(
            "sed 1,1d | tr '\t' '\n' | LC_ALL=C sort -g -k 1nr | uniq | head -n 1",
            input=mash_output,
            capture_output=True,
            text=True,
            shell=True,
        )

        max_divergence_str = divergence_values.stdout.strip()
        if not max_divergence_str:
            raise ValueError("No valid divergence found.")

        max_divergence = float(max_divergence_str)
        return math.ceil(100 - (max_divergence * 100)) - 2

    except subprocess.CalledProcessError as e:
        print(f"Error running Mash: {e}")
        sys.exit(1)


def prepare_fasta_and_analyze(
    input_dir: str, selected_files=None, output_yaml: str = "/output/params.yaml"
):
    """Prepares FASTA files for analysis, filters by selected files if provided."""
    merged_fasta_gz, downsampled_fasta, record_count, headers = rename_and_compress_fasta(
        input_dir, selected_files=selected_files
    )

    if not merged_fasta_gz:
        print("FASTA preparation failed.")
        sys.exit(1)

    percent_identity = run_mash_triangle(merged_fasta_gz)
    chosen_reference = random.choice(headers)

    results = {
        "path_to_input": merged_fasta_gz,
        "downsampled_fasta": downsampled_fasta,
        "number_of_haplotypes": record_count,
        "percent_identity": percent_identity,
        "reference_header": chosen_reference,
    }

    with open(output_yaml, "w") as yaml_file:
        yaml.dump(results, yaml_file, default_flow_style=False)

    print(f"Percentage identity found: {percent_identity}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare and analyse FASTA input for MineGraph."
    )
    parser.add_argument("input_dir",
                        help="Directory containing input FASTA files (mounted as /data)")
    parser.add_argument("output_dir", nargs="?", default="/output",
                        help="Output directory (mounted as /output); reserved for future use")
    parser.add_argument("--output_yaml", default="/output/params.yaml",
                        help="Path to write the output params YAML (default: /output/params.yaml)")
    parser.add_argument("--metadata", default=None,
                        help="Path to metadata CSV listing filenames to include (one per line)")
    parser.add_argument("selected_files", nargs="*",
                        help="Optional list of FASTA filenames to process (default: all)")

    args = parser.parse_args()

    # Explicit positional files take priority; fall back to metadata CSV; then use all files.
    selected: list | None = args.selected_files or None
    if not selected and args.metadata:
        try:
            with open(args.metadata) as fh:
                selected = [line.strip() for line in fh if line.strip()]
            if not selected:
                print(f"[WARN] --metadata file '{args.metadata}' is empty; using all files.")
                selected = None
            else:
                print(f"[INFO] Loaded {len(selected)} filenames from metadata: {args.metadata}")
        except OSError as exc:
            print(f"[ERROR] Cannot read --metadata file: {exc}", file=sys.stderr)
            sys.exit(1)

    prepare_fasta_and_analyze(
        args.input_dir,
        selected_files=selected,
        output_yaml=args.output_yaml,
    )
