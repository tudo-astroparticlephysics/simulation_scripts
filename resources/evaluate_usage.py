import re
import glob
import numpy as np
import warnings
import click
import os


# Ignore RuntimeWarning
warnings.filterwarnings("ignore", category=RuntimeWarning)


def extract_values_from_output_logs(log_file_path):
    with open(log_file_path, "r") as file:
        log_content = file.read()

    # Define the pattern to search for "That took " followed by a number
    pattern = re.compile(r"That took (\d+(\.\d+)?)")

    # Use findall to extract all matching values
    matches = pattern.findall(log_content)

    if len(matches) == 1:
        return matches[0][0]
    elif len(matches) < 1:
        return None
    else:
        raise ValueError(f"Found multiple matches in {log_file_path}")


def extract_values_from_job_logs(log_file_path, only_one_line=False):
    with open(log_file_path, "r") as file:
        log_content = file.read()

    # Define the pattern to search for
    mem_line = re.findall(r"\t   Memory.*", log_content, re.MULTILINE)
    disk_line = re.findall(r"\t   Disk.*", log_content, re.MULTILINE)
    num_pattern = re.compile("(\d+(\.\d+)?)")

    if len(mem_line) > 1 and only_one_line:
        raise ValueError(f"Found {len(mem_line)} memory lines in {log_file_path}")

    if len(disk_line) > 1 and only_one_line:
        raise ValueError(f"Found {len(disk_line)} disk lines in {log_file_path}")

    if len(mem_line) == 0:
        mem_usage = None
    else:
        matches_mem = num_pattern.findall(mem_line[-1])
        if len(matches_mem) == 3:
            mem_usage = matches_mem[0][0]
        else:
            mem_usage = None

    if len(disk_line) == 0:
        disk_usage = None
    else:
        matches_disk = num_pattern.findall(disk_line[-1])
        if len(matches_disk) == 3:
            disk_usage = matches_disk[0][0]
        else:
            disk_usage = None

    return mem_usage, disk_usage


def process_all_out_files(directory_path):
    # Use glob to get a list of all .out files in the specified directory
    out_files = glob.glob(f"{directory_path}/*.out")
    log_files = glob.glob(f"{directory_path}/*.log")

    # Process each .out file
    times = []
    for out_file in out_files:
        extracted_values = extract_values_from_output_logs(out_file)
        # print(f"File: {out_file}, Time: {extracted_values}")
        if extracted_values is not None:
            times.append(float(extracted_values))

    disk_usage = []
    mem_usage = []
    for log_file in log_files:
        mem, disk = extract_values_from_job_logs(log_file)
        if mem is not None:
            mem_usage.append(float(mem))
        if disk is not None:
            disk_usage.append(float(disk))

    # convert from KB to MB
    disk_usage = np.array(disk_usage) / 1024

    for name, unit, values in zip(
        ["Time", "Disk", "Memory"], ["s", "MB", "MB"], [times, disk_usage, mem_usage]
    ):
        if len(values) == 0:
            print(f"    {name:>12s}: No values found.")
            continue
        print(
            f"    {name:>12s}: {np.mean(values):6.0f} ± {np.std(values):<6.0f} {unit:<2s}"
            f" [{np.min(values):6.0f} {unit:>2s}, {np.max(values):>6.0f} {unit:>2s}]"
            f" (N: {len(values)})"
        )


def get_file_size(file_path):
    size_in_bytes = os.path.getsize(file_path)
    size_in_kb = size_in_bytes / 1024.0
    return size_in_kb


def check_disk_usage(directory_path, pattern="*/*.i3*"):
    files = glob.glob(os.path.join(directory_path, pattern))

    file_sizes = []
    for file_name in files:
        file_path = os.path.join(directory_path, file_name)
        if os.path.isfile(file_path):
            file_size = get_file_size(file_path)
            # print(f"File: {file_name}, Size: {file_size/1e3:.2f} MB")
            file_sizes.append(file_size)
    print(
        f"    {'File size':>12s}: {np.mean(file_sizes)/1e3:6.2f} ± "
        f"{np.std(file_sizes)/1e3:6.2f} MB (#files: {len(file_sizes)})"
    )


@click.command()
@click.argument("scratch_dir", type=click.Path(exists=True))
@click.option(
    "--out_folder",
    "-o",
    default=None,
    help="Directory where the output files are stored..",
)
@click.option(
    "--steps",
    "-s",
    multiple=True,
    default=[],
    help="The simulation step numbers to be analyzed. If None, all steps are analyzed.",
)
@click.option(
    "--datasets",
    "-d",
    multiple=True,
    default=[],
    help="The datasets to be analyzed. If None, all datasets are analyzed.",
)
@click.option(
    "--pattern",
    "-p",
    default="{dataset}_step_{step}_*",
    help="The pattern to be used to find the directories. Default: {dataset}_step_{step}_*",
)
def main(scratch_dir, out_folder, steps, datasets, pattern):
    # find datasets
    if len(datasets) == 0:
        dataset_candidates = [d.split("_")[0] for d in os.listdir(scratch_dir)]
        datasets = sorted(set([int(d) for d in dataset_candidates if d.isdigit()]))

    # find steps
    if len(steps) == 0:
        step_candidates = [d.split("_")[2] for d in os.listdir(scratch_dir)]
        steps = sorted(set([int(d) for d in step_candidates if d.isdigit()]))

    for dataset in datasets:
        print(f"Dataset: --- {dataset} ---")
        for step in steps:
            # find directory for this step and dataset
            pattern_i = pattern.format(dataset=dataset, step=step)
            directory = glob.glob(os.path.join(scratch_dir, pattern_i))
            assert (
                len(directory) == 1
            ), f"Found {len(directory)} directories for dataset {dataset} and step {step}."
            directory = directory[0]

            step_name = "_".join(os.path.basename(directory).split("_")[3:])

            print(f"  Step: {step} [{step_name}]")
            process_all_out_files(os.path.join(directory, "logs"))

            if out_folder is not None:
                data_dir = os.path.join(out_folder, f"{dataset}/step_{step}_*")
                check_disk_usage(data_dir)


if __name__ == "__main__":
    main()
