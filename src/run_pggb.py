import sys

import yaml
import subprocess
import os

def run_pggb(yaml_file="/output/params.yaml", threads=16, resume=False,
             transclose_batch=None, seqwish_temp_dir=None):
    # Load parameters from the YAML file
    try:
        with open(yaml_file, "r") as f:
            params = yaml.safe_load(f)

        # Extract parameters with defaults if not present
        path_to_input = params.get("path_to_input")
        number_of_haplotypes = params.get("number_of_haplotypes", 9)
        percent_identity = params.get("percent_identity", 90)
        segment_length = params.get("segment_length", 5000)
        vcf_reference = params.get("reference_header", "")

        # Convert segment_length to string with "k" if needed (e.g., 5000 -> "5k")
        if segment_length >= 1000:
            segment_length = f"{int(segment_length / 1000)}k"
        else:
            segment_length = str(segment_length)

        # Ensure the output directory exists
        output_dir = "/output/pggb_output"
        os.makedirs(output_dir, exist_ok=True)

        # Construct the pggb command
        pggb_command = [
            "pggb",
            "-i", path_to_input,
            "-o", output_dir,
            "-n", str(number_of_haplotypes),
            "-t", str(threads),
            "-p", str(percent_identity),
            "-s", segment_length,
            "-V", vcf_reference,
            "-m", "-M",
        ]

        # -B: transitive closure batch size in bp (default 10000000).
        # Increase to match your total input size so all sorting fits in one pass.
        # e.g. 709 chloroplast genomes ~150kb each = ~106Mbp → use 500000000 (500Mbp).
        if transclose_batch is not None:
            pggb_command += ["-B", str(transclose_batch)]

        # -D: temp dir for seqwish's disk-backed mmap sort files.
        # On AWS, point this at a tmpfs or NVMe mount, NOT the EBS output volume.
        if seqwish_temp_dir is not None:
            pggb_command += ["-D", seqwish_temp_dir]

        if resume:
            pggb_command.append("--resume")

        # Run the pggb command
        print("Running pggb with command:", " ".join(pggb_command))
        subprocess.run(pggb_command, check=True)
        print(f"pggb run completed. Output saved to {output_dir}")

    except FileNotFoundError:
        print(f"Error: {yaml_file} not found.")
    except yaml.YAMLError:
        print("Error reading the YAML file.")
    except subprocess.CalledProcessError:
        print("Error running pggb.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    import argparse as _ap
    _p = _ap.ArgumentParser(add_help=False)
    _p.add_argument("threads")
    _p.add_argument("--resume", action="store_true", default=False)
    _p.add_argument("--transclose-batch", dest="transclose_batch", type=int, default=None,
                    help="seqwish transitive closure batch size in base pairs (pggb -B). "
                         "Default: 10000000. Increase to reduce disk sort passes on high-RAM machines. "
                         "Example: 500000000 for 709 chloroplast genomes.")
    _p.add_argument("--seqwish-temp-dir", dest="seqwish_temp_dir", default=None,
                    help="Directory for seqwish mmap temp files (pggb -D). "
                         "Point to /tmp (tmpfs) or a local NVMe mount to avoid EBS I/O. "
                         "Example: /tmp")
    _a, _ = _p.parse_known_args()
    run_pggb(threads=int(_a.threads), resume=_a.resume,
             transclose_batch=_a.transclose_batch,
             seqwish_temp_dir=_a.seqwish_temp_dir)
