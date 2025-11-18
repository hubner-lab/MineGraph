import os
import sys
import pandas as pd
import subprocess
import argparse
import time
import resource
import platform
import os
import glob

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


def run_workflow(data_dir, output_dir, metadata_file, threads=16, tree_pars=10, tree_bs=10, quantile=0.25, top_n=50,
                 sq_view=False, w_size=1000, mode="all", phyl_tree=True, plots=True, only_stats=False,tree_type='graph',
                 softcore_threshold=90):
    try:
        machine = platform.machine()
        if machine == "x86_64":
            docker_platform = "linux/amd64"
        elif machine in ["arm64", "aarch64"]:
            docker_platform = "linux/arm64"
        else:
            docker_platform = None

        if not os.path.exists(metadata_file):
            print(f"[ERROR] Metadata file {metadata_file} not found.")
            sys.exit(1)

        # Read metadata file (CSV or XLSX)
        if metadata_file.endswith(".csv"):
            metadata = pd.read_csv(metadata_file, header=None)
        elif metadata_file.endswith(".xlsx"):
            metadata = pd.read_excel(metadata_file, header=None)
        else:
            print("[ERROR] Metadata file must be .csv or .xlsx")
            sys.exit(1)

        print(f"[INFO] Loaded metadata with {metadata.shape[0]} rows and {metadata.shape[1]} columns.")

        os.makedirs(output_dir, exist_ok=True)
        fasta_files = metadata.iloc[:, 0].tolist()
        if not fasta_files:
            print("[ERROR] No FASTA files found in the data directory.")
            sys.exit(1)
        
        # Step 1: Prepare FASTA input
        if mode in ["all", "extract-tr", "construct-graph"]:
            print("[STEP 1] Preparing FASTA input...")
            prepare_command = [
                                  "docker", "run", "--rm","--platform", f"{docker_platform}", "-v", f"{os.path.abspath(data_dir)}:/data",
                                  "-v", f"{os.path.abspath(output_dir)}:/output",
                                  "rakanhaib/opggb", "python", "/src/prepare_and_mash_input.py", "/data"
                              ] + fasta_files
            subprocess.run(prepare_command, check=True)

        # Step 2: Run RepeatMasker
        if mode in ["all", "extract-tr", "construct-graph"]:
            print("[STEP 2] Running RepeatMasker...")
            repeatmask_command = [
                "docker", "run", "--rm","--platform", f"{docker_platform}", "-v", f"{os.path.abspath(output_dir)}:/data",
                "pegi3s/repeat_masker", "bash", "-c",
                f"RepeatMasker /data/downsampled_panSN_output.fasta -no_is -pa {threads} -s -gff"
            ]
            subprocess.run(repeatmask_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Step 3: Extract longest TR and update parameters
        if mode in ["all", "extract-tr", "construct-graph"]:
            print("[STEP 3] Extracting longest TR and updating parameters...")
            extract_command = [
                "docker", "run", "--rm","--platform", f"{docker_platform}", "-v", f"{os.path.abspath(output_dir)}:/data",
                "rakanhaib/opggb", "python", "/src/run_repeatmask.py"
            ]
            subprocess.run(extract_command, check=True)

        if mode == "extract-tr":
            print("[INFO] Extract-only mode complete.")
            return

        # Step 4: Run PGGB
        print(f"[STEP 4] Running PGGB with {threads} threads...")
        pggb_command = [
            "docker", "run", "--rm","--platform", f"{docker_platform}", "-v", f"{os.path.abspath(output_dir)}:/output",
            "rakanhaib/opggb", "python", "/src/run_pggb.py", str(threads)
        ]
        subprocess.run(pggb_command, check=True)

        if mode == "construct-graph":
            print("[INFO] Graph construction complete. Skipping statistics.")
            return

        
        # Step 5: Run Graph Statistics Analysis
        print("[STEP 5] Performing graph statistical analysis...")
        stats_command = [
            "docker", "run", "--rm","--platform", f"{docker_platform}",
            "-v", f"{os.path.abspath(output_dir)}:/data",
            "rakanhaib/opggb",
            "python", "/src/run_stats.py",
            "--threads", str(threads),
            "--tree_pars", str(tree_pars),
            "--tree_bs", str(tree_bs),
            "--input_dir", "/data/pggb_output",
            "--output_dir", "/data",
            "--quantile", str(quantile / 100),
            "--top_n", str(top_n),
            "--window_size", str(w_size),
            "--tree_type", str(tree_type),
            "--softcore_threshold", str(softcore_threshold),
        ]

        # Add new statistical arguments
        if phyl_tree:
            stats_command.append("--phyl-tree")
        if plots:
            stats_command.append("--plots")
        if only_stats:
            stats_command.append("--only-stats")

        subprocess.run(stats_command, check=True)

        if sq_view:
            stp_command = [
                "docker", "run","--platform", f"{docker_platform}", "-it", "--rm",
                "-v", f"{os.path.abspath(output_dir)}/MineGraph_output/:/data",
                "-p", "3210:3000",
                "rakanhaib/sequencetubemap:latest",
                "bash", "-c",
                "./scripts/prepare_vg.sh /data/gfa_to_vg.vg && npm start serve"
            ]
            subprocess.run(stp_command, check=True)
        print("[INFO] Statistical analysis complete.")
        end_time = time.time()
        usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)
        elapsed_time = end_time - start_time
        max_memory_usage_mb = 0
        print("\n===== PIPELINE RESOURCE USAGE =====")
        print(f"Elapsed wall clock time: {elapsed_time:.2f} seconds")
        # Define the log directory path
        log_dir = os.path.join(os.path.abspath(os.path.expanduser(args.output_dir)), "pggb_output")
        log_files = glob.glob(os.path.join(log_dir, "*.log"))
        total_cpu = 0.0
        rax_log_files = glob.glob(os.path.join(log_dir, "*.log"))
        try:
            with open(log_files[0], 'r') as f:
                for line in f:
                    if "max memory" in line:
                        # Extract the time value (7th column), remove 's' and convert to float
                        parts = line.split()
                        if len(parts) >= 7:
                            time_str = parts[6].replace('s', '')
                            mem_str = parts[8].replace('Kb', '')
                            try:
                                total_cpu += float(time_str)
                                if int(mem_str) > max_memory_usage_mb:
                                    max_memory_usage_mb = int(mem_str)
                            except ValueError:
                                print(f"Could not convert time value in {log_file}: {parts[6]}")
            if phyl_tree:
                rax_log_dir = os.path.join(os.path.abspath(os.path.expanduser(args.output_dir)), 'graph.raxml.log')
                with open(rax_log_dir, 'r') as rf:
                    try:
                        rf_lines = rf.readlines()
                        total_cpu += float(rf_lines[-2].split()[2])
                    except ValueError:
                        print(f"Could not convert mem value in {rax_log_files}")
        except IOError as e:
            print(f"Error reading {log_file}: {e}")
        print(f"Total CPU time (seconds): {total_cpu:.2f}")
        print(f"Max memory usage: {max_memory_usage_mb/1024:.2f} MB")
        print("====================================\n")
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Exiting...")
        sys.exit(-1)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}: {e}")
        sys.exit(e.returncode)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MineGraph Workflow with different processing options.")
    parser.add_argument("--data_dir", required=True, help="Directory containing input FASTA files.")
    parser.add_argument("--output_dir", required=True, help="Directory for saving outputs.")
    parser.add_argument("--metadata", required=True,
                        help="Metadata file (CSV or XLSX) containing additional information.")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use (default: 16).")
    parser.add_argument("--tree_pars", type=int, default=10, help="Number of parsimonious trees (default: 10).")
    parser.add_argument("--tree_bs", type=int, default=10, help="Number of bootstrap trees (default: 10).")
    parser.add_argument("--quantile", type=float, default=25, help="Consensus node percentage threshold (default: 25).")
    parser.add_argument("--top_n", type=int, default=50, help="Top N nodes to visualize (default: 50).")
    parser.add_argument("--window_size", type=int, default=1000, help="Sliding window size (default: 1000).")
    parser.add_argument("--sq_view", action="store_true", help="Run SequenceTubeMap viewer.")
    parser.add_argument("--mode", choices=["all", "extract-tr", "construct-graph"], default="all",
                        help="Select which steps to run.")
    parser.add_argument("--phyl-tree", action="store_true", help="Generate an MSA and a phylogenetic tree.")
    parser.add_argument("--plots", action="store_true", help="Generate all statistical plots.")
    parser.add_argument("--only-stats", action="store_false", help="Generate only the statistical files (xlsx).")
    parser.add_argument(
        "--tree_type",
        choices=["graph", "msa"],
        default="graph",
        help="Which phylogenetic tree to build if --phyl_tree is set: 'graph' (default, from ODGI matrix) or 'msa' (RAxML on MSA)."
    )
    parser.add_argument(
        "--sc_th",
        type=int,
        default=5,
        help="Percentage threshold defining the cutoff for unique and core nodes (default: 5%)."
             " Nodes present in ≤ threshold% of haplotypes are 'unique', and nodes in ≥ (100 - threshold)% are 'core'."
    )

    args = parser.parse_args()
    start_time = time.time()
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)

    run_workflow(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        metadata_file=args.metadata,
        threads=args.threads,
        tree_pars=args.tree_pars,
        tree_bs=args.tree_bs,
        quantile=args.quantile,
        top_n=args.top_n,
        w_size=args.window_size,
        sq_view=args.sq_view,
        mode=args.mode,
        phyl_tree=args.phyl_tree,
        plots=args.plots,
        only_stats=args.only_stats,
        tree_type=args.tree_type,
        softcore_threshold=args.sc_th
    )

