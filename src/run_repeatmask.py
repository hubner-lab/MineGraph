import yaml


def extract_longest_single_repeat(out_file):
    """Returns the length of the longest individual repeat."""
    longest_length = 0

    with open(out_file, "r") as f:
        for line in f:
            if line.startswith(("SW", "score")) or not line.strip():
                continue

            cols = line.split()
            if len(cols) >= 15:  # Match your original column requirement
                try:
                    # Using columns[11] and [12] to match your original logic ("position in repeat")
                    begin, end = int(cols[11]), int(cols[12])
                    repeat_length = end - begin + 1
                    if repeat_length > longest_length:
                        longest_length = repeat_length
                except (ValueError, IndexError):
                    continue

    return longest_length


def extract_longest_merged_region(out_file):
    """Returns the length of the longest merged repeat region."""
    intervals = []

    with open(out_file, "r") as f:
        for line in f:
            if line.startswith(("SW", "score")) or not line.strip():
                continue

            cols = line.split()
            if len(cols) >= 15:
                try:
                    # Using query sequence coordinates (cols[5] and [6]) for merging
                    start, end = int(cols[5]), int(cols[6])
                    intervals.append((start, end))
                except (ValueError, IndexError):
                    continue

    if not intervals:
        return 0

    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged:
            merged.append([start, end])
        else:
            last_start, last_end = merged[-1]
            if start <= last_end:
                merged[-1][1] = max(last_end, end)
            else:
                merged.append([start, end])

    return max(end - start + 1 for start, end in merged)


def update_yaml_with_segment_length(yaml_file, segment_length):
    """Updates the YAML file with the calculated segment length."""
    try:
        with open(yaml_file, "r") as yf:
            data = yaml.safe_load(yf) or {}

        data["segment_length"] = segment_length

        with open(yaml_file, "w") as yf:
            yaml.dump(data, yf, default_flow_style=False)

        print(f"Updated segment_length in {yaml_file}: {segment_length}")
    except Exception as e:
        print(f"Error updating YAML file: {e}")


def extract_longest_repeat(out_file, yaml_file="/data/params.yaml"):
    """Wrapper function matching your original interface."""
    try:
        # Get both lengths
        single_length = extract_longest_single_repeat(out_file)
        merged_length = extract_longest_merged_region(out_file)

        # Use max of both (or single_length if you prefer)
        longest_repeat_length = max(single_length, merged_length)

        # Apply your 1.2x multiplier
        segment_length = int(longest_repeat_length * 1.2)

        # Update YAML
        update_yaml_with_segment_length(yaml_file, segment_length)

    except Exception as e:
        print(f"Error processing {out_file}: {e}")


if __name__ == "__main__":
    # Maintain your original usage
    extract_longest_repeat("/data/downsampled_panSN_output.fasta.out")