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
        python run_workflow.py --data_dir <data_dir> --output_dir <output_dir> [OPTIONS]

    REQUIRED ARGUMENTS:
        --data_dir    : Directory containing input FASTA files.
        --output_dir  : Directory where results will be saved.

    OPTIONAL ARGUMENTS:
        --metadata    : File listing FASTA files to process (CSV/XLSX). 
                        If not provided, all FASTA files in the data directory will be processed.
        --threads     : Number of threads to use for PGGB and other parallel tasks (default: 16).
        --tree_pars   : Number of parsimonious trees to generate (default: 10).
        --tree_bs     : Number of bootstrap trees to generate (default: 10).
        --quantile    : Consensus nodes percentage of presence (default: 50). 
                        Example: 100 means the nodes appear in 100% of the paths.
        --top_n       : Number of top nodes (by size) to visualize (default: 1000).

    EXAMPLES:
        Example 1: Run with selected FASTA files and 64 threads
            $ python run_workflow.py --data_dir /path/to/data --output_dir /path/to/output --metadata selected_files.csv --threads 64

    LET'S GET STARTED!
    """)


def run_workflow(data_dir, output_dir, metadata=None, threads=16, tree_pars=10, tree_bs=10, quantile=0.25, top_n=50, sq_view=False,w_size=1000):
    """
    Run the full workflow with Docker, starting from the given directory.
    Accepts either a file list or all files in data_dir if no file is provided.

    Args:
        data_dir (str): Directory with FASTA files.
        output_dir (str): Directory where results will be saved.
        metadata (str, optional): CSV or XLSX file with list of FASTA files to process.
        threads (int): Number of threads to use for PGGB and other parallel tasks.
        :param w_size:
        :param sq_view:
        :param top_n:
        :param threads:
        :param metadata:
        :param output_dir:
        :param data_dir:
        :param quantile:
        :param tree_bs: an argument used to be passed for Raxmel bootstraps --bs-trees
        :param tree_pars: an argument used to be passed for Raxmel parsimonious trees --tree
    """
    try:

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
        """
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
        subprocess.run(repeatmask_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
        """
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
            "--top_n", "{}".format(top_n),
            "--window_size", "{}".format(w_size)
        ]
        subprocess.run(stats_command, check=True)
        if sq_view:
            stp_command = [
                "docker", "run", "-it", "--rm",
                "-v", f"{os.path.abspath(output_dir)}/MineGraph_output/:/data",
                "-p", "3210:3000",
                "rakanhaib/sequencetubemap:latest",
                "bash", "-c",
                "./scripts/prepare_vg.sh /data/gfa_to_vg.vg && npm start serve"
            ]
            subprocess.run(stp_command, check=True)
        print("[INFO] Statistical analysis completed.")

        print(
            "\n🎉 [WORKFLOW COMPLETE] All steps finished successfully. Results are in the specified output directory. 🎉")
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Exiting...")
        sys.exit(-1)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}: {e}")
        sys.exit(e.returncode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MineGraph Workflow: Process input FASTA files and visualize sequence graphs."
    )

    parser.add_argument("--data_dir", required=True, help="Directory containing input FASTA files.")
    parser.add_argument("--output_dir", required=True, help="Directory for saving outputs.")
    parser.add_argument("--metadata", required=True, help="File listing FASTA files to process (CSV/XLSX).")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use (default: 16).")
    parser.add_argument("--tree_pars", type=int, default=10,
                        help="Number of parsimonious trees to generate (default: 10).")
    parser.add_argument("--tree_bs", type=int, default=10, help="Number of bootstrap trees (default: 10).")
    parser.add_argument("--quantile", type=int, default=50,
                        help="Consensus nodes percentage of presence. "
                             "100 means the nodes appeared in all paths (default: 50).")
    parser.add_argument("--top_n", type=int, default=1000,
                        help="Number of top nodes by size to visualize (default: 1000).")
    parser.add_argument("--view", action="store_true",
                        help="Launch SequenceTubeMap if this flag is provided.")
    parser.add_argument("--window_size", type=int, default=100,
                        help="Sliding window size for analysis (default: 100).")

    args = parser.parse_args()

    run_workflow(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        metadata=args.metadata,
        threads=args.threads,
        tree_pars=args.tree_pars,
        tree_bs=args.tree_bs,
        quantile=args.quantile / 100,
        top_n=args.top_n,
        sq_view=args.view,
        w_size=args.window_size
    )
