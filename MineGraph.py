import os
import sys
import pandas as pd
import subprocess
import argparse


def print_help():
    print(r"""
    ============================================================
    || WELCOME TO THE MineGraph - GRAPH ANALYSIS WORKFLOW      ||
    ============================================================
    |                                                          |
    |     This workflow will take you through:                 |
    |     - FASTA preparation and compression                  |
    |     - Repeat Masker analysis for Tandem Repeats          |
    |     - PGGB alignment and graph generation                |
    |     - Graph statistical analysis                         |
    |                                                          |
    |                   Developed witH ❤️                      |
    ============================================================
    """)

    print("""
    HOW TO USE THIS WORKFLOW:

    Usage:
        python run_workflow.py --data_dir <data_dir> --output_dir <output_dir> [--metadata <metadata>] [--threads <threads>]

    Arguments:
        --data_dir    : Directory containing the raw FASTA files.
        --output_dir  : Directory where results will be saved.
        --metadata    : (Optional) CSV or XLSX file with list of FASTA files to process.
                        If not provided, all FASTA files in the <data_dir> will be used.
                        File must have one column with FASTA filenames.
        --threads     : (Optional) Number of threads to use for PGGB and other parallel tasks.
                        Default is 16.

    EXAMPLES:
        Example 1: Run with selected FASTA files with 64 threads
            $ python run_workflow.py --data_dir /path/to/data --output_dir /path/to/output --metadata selected_files.csv --threads 64

        Example 2: Run with all FASTA files in the data directory
            $ python run_workflow.py --data_dir /path/to/data --output_dir /path/to/output 

    LET'S GET STARTED!
    """)


def run_workflow(data_dir, output_dir, metadata=None, threads=16, tree_pars=10, tree_bs=10, quantile=0.25, top_n=50):
    """
    Run the full workflow with Docker, starting from the given directory.
    Accepts either a file list or all files in data_dir if no file is provided.

    Args:
        data_dir (str): Directory with FASTA files.
        output_dir (str): Directory where results will be saved.
        metadata (str, optional): CSV or XLSX file with list of FASTA files to process.
        threads (int): Number of threads to use for PGGB and other parallel tasks.
        :param threads:
        :param metadata:
        :param output_dir:
        :param data_dir:
        :param quantile:
        :param tree_bs: an argument used to be passed for Raxmel bootstraps --bs-trees
        :param tree_pars: an argument used to be passed for Raxmel parsimonious trees --tree
    """
    os.makedirs(output_dir, exist_ok=True)

    # Determine selected files based on metadata or process all files in data_dir
    if metadata:
        if metadata.endswith(".csv"):
            fasta_files_df = pd.read_csv(metadata, header=None)
        elif metadata.endswith(".xlsx"):
            fasta_files_df = pd.read_excel(metadata, header=None)
        else:
            print("[ERROR] Unsupported file format. Please use CSV or XLSX.")
            sys.exit(1)

        if fasta_files_df.shape[1] != 1:
            print("[ERROR] Metadata file must contain only one column with FASTA file names.")
            sys.exit(1)

        fasta_files = fasta_files_df.iloc[:, 0].tolist()
        print(f"[INFO] Processing {len(fasta_files)} FASTA files from {metadata}")
    else:
        fasta_files = [f for f in os.listdir(data_dir) if f.endswith(".fasta")]
        print(f"[INFO] Processing all FASTA files in {data_dir}.")

    # Step 1: Prepare FASTA input
    prepare_command = [
                          "docker", "run", "--rm", "-v", f"{os.path.abspath(data_dir)}:/data",
                          "-v", f"{os.path.abspath(output_dir)}:/output",
                          "rakanhaib/opggb", "python", "/prepare_and_mash_input.py", "/data"
                      ] + fasta_files
    subprocess.run(prepare_command, check=True)

    # Step 2: Run RepeatMasker on the downsampled FASTA file
    print("[STEP 2/5] Running RepeatMasker on downsampled FASTA...")
    repeatmask_command = [
    "docker", "run", "--rm", "-v", f"{os.path.abspath(output_dir)}:/data",
    "pegi3s/repeat_masker", "bash", "-c",
    f"RepeatMasker /data/downsampled_panSN_output.fasta -pa {threads} -no_is -s"]
    subprocess.run(repeatmask_command, check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("[INFO] RepeatMasker analysis completed.")

    # Step 3: Run run_repeatmask.py inside Docker to extract the longest TR and update params.yaml
    print("[STEP 3/5] Extracting longest tandem repeat and updating parameters...")
    extract_command = [
        "docker", "run", "--rm", "-v", f"{os.path.abspath(output_dir)}:/data",
        "rakanhaib/opggb", "python", "/run_repeatmask.py"
    ]
    subprocess.run(extract_command, check=True)
    print("[INFO] Longest tandem repeat extraction completed and params.yaml updated.")

    # Step 4: Run run_pggb.py inside Docker using the specified number of threads
    print(f"[STEP 4/5] Running PGGB with {threads} threads...")
    pggb_command = [
        "docker", "run", "--rm", "-v", f"{os.path.abspath(output_dir)}:/output",
        "rakanhaib/opggb", "python", "/run_pggb.py", str(threads)
    ]
    subprocess.run(pggb_command, check=True)
    print("[INFO] PGGB alignment and graph generation completed.")

    # Step 5: Run run_stats.py inside Docker for statistical analysis
    print("[STEP 5/5] Performing statistical analysis on generated graph and alignments...")
    stats_command = [
        "docker", "run", "--rm",
        "-v", f"{os.path.abspath(output_dir)}:/data",
        "rakanhaib/opggb",
        "python", "/run_stats.py",
        "--threads", "{}".format(threads),
        "--tree_pars", "{}".format(tree_pars),
        "--tree_bs", "{}".format(tree_bs),
        "--input_dir", "/data/pggb_output",
        "--output_dir", "/data/MineGraph_output",
        "--quantile", "{}".format(quantile),
        "--top_n", "{}".format(top_n)

    ]
    subprocess.run(stats_command, check=True)
    print("[INFO] Statistical analysis completed.")

    print("\n🎉 [WORKFLOW COMPLETE] All steps finished successfully. Results are in the specified output directory. 🎉")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MineGraph Workflow")
    parser.add_argument("--data_dir", required=True, help="Directory containing input FASTA files")
    parser.add_argument("--output_dir", required=True, help="Directory for saving outputs")
    parser.add_argument("--metadata", required=True, help="File listing FASTA files to process (CSV/XLSX)")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads (default: 16)")
    parser.add_argument("--tree_pars", type=int, default=10, help="Number of parsimonious trees (default: 10)")
    parser.add_argument("--tree_bs", type=int, default=10, help="Number of bootstrap trees (default: 10)")
    parser.add_argument("--quantile", type=int, default=50, help="Consensus nodes percentage of presence,"
                                                                 " 100 means the nodes appeared in 100% of the paths, default is 50")
    parser.add_argument("--top_n", type=int, default="1000", help="top N nodes sizes to be visualized, default is 1000")

    args = parser.parse_args()

    run_workflow(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        metadata=args.metadata,
        threads=args.threads,
        tree_pars=args.tree_pars,
        tree_bs=args.tree_bs,
        quantile=args.quantile/100,
        top_n=args.top_n
    )
