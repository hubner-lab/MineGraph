import sys

import yaml
import subprocess
import os

def run_pggb(yaml_file="/output/params.yaml", threads=16):
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
            "-m", "-M"

        ]

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
    threads = int(sys.argv[1])
    run_pggb(threads=threads)
