#!/usr/bin/env python3
"""
Author : conrad <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-09-09
Purpose: Create Trimmomatic jobs for ARC.
"""

import argparse
import os.path
import sys


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Create Trimmomatic or FastQC jobs for ARC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # ----------------------------------------------
    # General options
    # ----------------------------------------------
    parser.add_argument(
        "-an",
        "--analysis",
        help="Analysis to run, one of: trimmomatic or fastqc. Required for FastQC and Trimmomatic.",
        type=str,
        metavar="str",
        required=True,
    )

    parser.add_argument(
        "-nt",
        "--num_threads",
        help="Number of threads to use for job. Set based on partitions desired. Default=8 appropriate for partitions single,lattice.",
        type=int,
        default=8,
        metavar="int",
    )

    parser.add_argument(
        "-mx",
        "--memmax",
        help="Max memory to use for job in MB; default=all available (0). Required for FastQC and Trimmomatic.",
        type=int,
        default=0,
        metavar="int",
    )

    parser.add_argument(
        "-mt",
        "--maxtime",
        help="max time for job to run. Required for FastQC and Trimmomatic.",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-p",
        "--partitions",
        help='ARC partitions to use, comma-separated list; default="single,lattice". Required for FastQC and Trimmomatic.',
        type=str,
        default="single,lattice",
        metavar="str",
    )

    parser.add_argument(
        "-i",
        "--isolate_list",
        help="path to isolate list on local computer; one isolate per line in list. Required for FastQC and Trimmomatic.",
        metavar="file",
        required=True,
    )

    parser.add_argument(
        "-cs",
        "--chunk_size",
        help="number of isolates (chunk) to analyze in one job file. Required for FastQC and Trimmomatic. Recommended size for Unicycler = 1.",
        type=int,
        metavar="int",
        required=True,
    )

    parser.add_argument(
        "-fe",
        "--fastq_ending",
        help="fastq file ending, one of: .fastq, .fq, .fastq.gz, or .fq.gz. Required for FastQC and Trimmomatic.",
        type=str,
        default=".fastq.gz",
        metavar="str",
    )

    parser.add_argument(
        "-d",
        "--input_dir",
        help="input directory on ARC containing fastq files to be trimmed. Default = current dir. Required for FastQC and Trimmomatic.",
        type=str,
        metavar="path/to/dir",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="output directory on ARC where trimmed files will be placed. Default = current dir. Required for FastQC and Trimmomatic.",
        type=str,
        metavar="path/to/dir",
        required=True,
    )

    parser.add_argument(
        "-j",
        "--jobs_dir",
        help="directory (local) to write job files to. default = current dir. Required for FastQC and Trimmomatic.",
        type=str,
        metavar="path/to/dir",
        required=True,
    )

    parser.add_argument(
        "-pre",
        "--prefix",
        help="job name prefix for job scripts, .out, and .err files. Required for FastQC and Trimmomatic.",
        type=str,
        metavar="str",
        required=True,
    )

    # ----------------------------------------------
    # Trimmomatic options
    # ----------------------------------------------

    parser.add_argument(
        "-a",
        "--adapters_file",
        help="path to adapters file on ARC. Required for Trimmomatic.",
        type=str,
        metavar="file",
    )

    parser.add_argument(
        "-ml",
        "--min_len",
        help="minimum length for trimmed reads, MINLEN argument in Trimmomatic. Required for Trimmomatic.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-r",
        "--read_len",
        help="read length to crop to, CROP argument in Trimmomatic. Required for Trimmomatic.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-q",
        "--min_qual",
        help="minimum quality in 4bp sliding window; SLIDINGWINDOW argument in Trimmomatic. Required for Trimmomatic.",
        type=int,
        metavar="int",
    )

    # ----------------------------------------------
    # Unicycler options
    # ----------------------------------------------

    parser.add_argument(
        "-mo",
        "--mode",
        help="Unicycler mode to run. Choose from bold, normal, or conservative. Required for Unicycler.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-dp",
        "--depth_filter",
        help="Unicycler depth filter to use. Recommended value = 0.25. Required for Unicycler.",
        type=float,
        metavar="float",
    )

    parser.add_argument(
        "-mc",
        "--min_contig_len",
        help="Minimum contig length Unicycler should produce. Required for Unicycler.",
        type=int,
        metavar="int",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Create Trimmomatic or FastQC jobs."""

    args = get_args()

    if args.analysis not in ["fastqc", "trimmomatic", "unicycler"]:
        print(
            "Error: analysis type specified is not valid. Analysis type (--analysis) must be either fastqc or trimmomatic."
        )
        sys.exit()

    if args.analysis == "trimmomatic":

        if args.adapters_file is None:
            print(
                "\nPlease provide a path to an adapters file for trimming using --adapters_file or -a. Required for creation of Trimmomatic jobs. Exiting.\n"
            )
            sys.exit()
        if args.min_len is None:
            print(
                "\nPlease provide a minimum read length to keep after trimming using --minlen or -ml. Required for creation of Trimmomatic jobs. Exiting.\n"
            )
            sys.exit()
        if args.read_len is None:
            print(
                "\nPlease provide an expected read length (i.e. read length from sequencing) using --read_length or -r. Required for creation of Trimmomatic jobs. Exiting.\n"
            )
            sys.exit()
        if args.min_qual is None:
            print(
                "\nPlease provide a minimum quality for trimming using --min_qual or -q. Required for creation of Trimmomatic jobs. Exiting.\n"
            )
            sys.exit()

        num_threads_arg = args.num_threads
        memmax_arg = args.memmax
        maxtime_arg = args.maxtime
        partitions_arg = args.partitions
        isolate_list_arg = args.isolate_list
        chunk_size_arg = args.chunk_size
        fastq_ending_arg = args.fastq_ending
        input_dir_arg = args.input_dir
        output_dir_arg = args.output_dir
        job_prefix_arg = args.prefix
        adapters_file_arg = args.adapters_file
        min_len_arg = args.min_len
        read_len_arg = args.read_len
        min_qual_arg = args.min_qual
        jobs_dir_arg = args.jobs_dir

        create_trimmomatic_jobs(
            num_threads_arg,
            memmax_arg,
            maxtime_arg,
            partitions_arg,
            isolate_list_arg,
            chunk_size_arg,
            fastq_ending_arg,
            input_dir_arg,
            output_dir_arg,
            job_prefix_arg,
            adapters_file_arg,
            min_len_arg,
            min_qual_arg,
            read_len_arg,
            jobs_dir_arg,
        )

    if args.analysis == "fastqc":
        num_threads_arg = args.num_threads
        memmax_arg = args.memmax
        maxtime_arg = args.maxtime
        partitions_arg = args.partitions
        isolate_list_arg = args.isolate_list
        chunk_size_arg = args.chunk_size
        fastq_ending_arg = args.fastq_ending
        input_dir_arg = args.input_dir
        output_dir_arg = args.output_dir
        job_prefix_arg = args.prefix
        jobs_dir_arg = args.jobs_dir

        create_fastqc_jobs(
            num_threads_arg,
            memmax_arg,
            maxtime_arg,
            partitions_arg,
            isolate_list_arg,
            chunk_size_arg,
            fastq_ending_arg,
            input_dir_arg,
            output_dir_arg,
            job_prefix_arg,
            jobs_dir_arg,
        )

    if args.analysis == "unicycler":

        if args.mode is None:
            print(
                "Please provide a mode for Unicycler using --mode. Required to create Unicycler jobs. Exiting."
            )
            sys.exit()
        if args.depth_filter is None:
            print(
                "Please provide a depth filter to Unicycler using --depth_filter. Required to create Unicycler jobs. Exiting."
            )
            sys.exit()
        if args.min_contig_len is None:
            print(
                "Please provide a minimum contig length using --min_contig_len. Required to create Unicycler jobs. Exiting."
            )
            sys.exit()

        num_threads_arg = args.num_threads
        memmax_arg = args.memmax
        maxtime_arg = args.maxtime
        partitions_arg = args.partitions
        isolate_list_arg = args.isolate_list
        chunk_size_arg = args.chunk_size
        fastq_ending_arg = args.fastq_ending
        input_dir_arg = args.input_dir
        output_dir_arg = args.output_dir
        job_prefix_arg = args.prefix
        jobs_dir_arg = args.jobs_dir
        mode_arg = args.mode
        depth_filter_arg = args.depth_filter
        min_contig_len_arg = args.min_contig_len

        create_unicycler_jobs(
            num_threads_arg,
            memmax_arg,
            maxtime_arg,
            partitions_arg,
            isolate_list_arg,
            chunk_size_arg,
            fastq_ending_arg,
            input_dir_arg,
            output_dir_arg,
            job_prefix_arg,
            jobs_dir_arg,
            mode_arg,
            depth_filter_arg,
            min_contig_len_arg,
        )


# --------------------------------------------------
def create_trimmomatic_jobs(
    num_threads_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    isolate_list_arg,
    chunk_size_arg,
    fastq_ending_arg,
    input_dir_arg,
    output_dir_arg,
    job_prefix_arg,
    adapters_file_arg,
    min_len_arg,
    min_qual_arg,
    read_len_arg,
    jobs_dir_arg,
):
    """Function to run Trimmomatic"""

    print(f'\nnum_threads_arg = "{num_threads_arg}"')
    print(f'memmax_arg = "{memmax_arg}"')
    print(f'maxtime_arg = "{maxtime_arg}"')
    print(f'partitions_arg = "{partitions_arg}"')
    print(f'isolate_list_arg = "{isolate_list_arg}"')
    print(f'chunk_size_arg = "{chunk_size_arg}"')
    print(f'fastq_ending_arg = "{fastq_ending_arg}"')
    print(f'input_dir_arg = "{input_dir_arg}"')
    print(f'output_dir_arg = "{output_dir_arg}"')
    print(f'job_prefix_arg = "{job_prefix_arg}"')
    print(f'adapters_file_arg = "{adapters_file_arg}"')
    print(f'min_len_arg = "{min_len_arg}"')
    print(f'read_len_arg = "{read_len_arg}"')
    print(f'min_qual_arg = "{min_qual_arg}"')
    print(f'jobs_dir_arg = "{jobs_dir_arg}"\n')

    # Read in isolate list:
    with open(isolate_list_arg, "r") as isolate_list_file:
        isolate_list = [isolate.strip() for isolate in isolate_list_file]

    # Create trimmomatic commands:
    trimmomatic_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(input_dir_arg, f"{isolate}_1{fastq_ending_arg}")
        input_r2 = os.path.join(input_dir_arg, f"{isolate}_2{fastq_ending_arg}")
        output_r1 = os.path.join(output_dir_arg, f"{isolate}_1{fastq_ending_arg}")
        output_r2 = os.path.join(output_dir_arg, f"{isolate}_2{fastq_ending_arg}")
        output_r1_unpaired = os.path.join(
            output_dir_arg, f"{isolate}_u_1{fastq_ending_arg}"
        )
        output_r2_unpaired = os.path.join(
            output_dir_arg, f"{isolate}_u_2{fastq_ending_arg}"
        )

        trimmomatic_cmd = (
            f"trimmomatic PE -threads {num_threads_arg} {input_r1} "
            f"{input_r2} {output_r1} {output_r1_unpaired} {output_r2} {output_r2_unpaired} "
            f"ILLUMINACLIP:{adapters_file_arg}:2:30:10:8:true CROP:{read_len_arg} "
            f"SLIDINGWINDOW:4:{min_qual_arg} MINLEN:{min_len_arg}\n"
        )

        trimmomatic_commands[isolate] = trimmomatic_cmd

    # Create chunked job files:
    chunks = [
        isolate_list[i : i + chunk_size_arg]
        for i in range(0, len(isolate_list), chunk_size_arg)
    ]

    for idx, chunk in enumerate(chunks):
        job_script = os.path.join(jobs_dir_arg, f"{idx}.{job_prefix_arg}.slurm")

        # Create job header:
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={maxtime_arg}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={num_threads_arg}\n"
            f"#SBATCH --partition={partitions_arg}\n"
            f"#SBATCH --mem={memmax_arg}\n"
            f"#SBATCH --output={idx}.{job_prefix_arg}.out\n"
            f"#SBATCH --error={idx}.{job_prefix_arg}.err\n"
            f"\nsource activate trimmomatic-env\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = trimmomatic_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


def create_fastqc_jobs(
    num_threads_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    isolate_list_arg,
    chunk_size_arg,
    fastq_ending_arg,
    input_dir_arg,
    output_dir_arg,
    job_prefix_arg,
    jobs_dir_arg,
):
    """Function to run fastqc"""

    print(f'\nnum_threads_arg = "{num_threads_arg}"')
    print(f'memmax_arg = "{memmax_arg}"')
    print(f'maxtime_arg = "{maxtime_arg}"')
    print(f'partitions_arg = "{partitions_arg}"')
    print(f'isolate_list_arg = "{isolate_list_arg}"')
    print(f'chunk_size_arg = "{chunk_size_arg}"')
    print(f'fastq_ending_arg = "{fastq_ending_arg}"')
    print(f'input_dir_arg = "{input_dir_arg}"')
    print(f'output_dir_arg = "{output_dir_arg}"')
    print(f'job_prefix_arg = "{job_prefix_arg}"')
    print(f'jobs_dir_arg = "{jobs_dir_arg}"')

    # Read in isolate list:
    with open(isolate_list_arg, "r") as isolate_list_file:
        isolate_list = [isolate.strip() for isolate in isolate_list_file]

    # Create fastqc commands:
    fastqc_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(input_dir_arg, f"{isolate}_1{fastq_ending_arg}")
        input_r2 = os.path.join(input_dir_arg, f"{isolate}_2{fastq_ending_arg}")

        fastqc_cmd_r1 = f"fastqc -t {num_threads_arg} -o {output_dir_arg} {input_r1}\n"
        fastqc_cmd_r2 = f"fastqc -t {num_threads_arg} -o {output_dir_arg} {input_r2}\n"

        fastqc_commands[isolate] = "\n".join([fastqc_cmd_r1, fastqc_cmd_r2])

    # Create chunked job files:
    chunks = [
        isolate_list[i : i + chunk_size_arg]
        for i in range(0, len(isolate_list), chunk_size_arg)
    ]

    for idx, chunk in enumerate(chunks):
        job_script = os.path.join(jobs_dir_arg, f"{idx}.{job_prefix_arg}.slurm")

        # Create job header:
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={maxtime_arg}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={num_threads_arg}\n"
            f"#SBATCH --partition={partitions_arg}\n"
            f"#SBATCH --mem={memmax_arg}\n"
            f"#SBATCH --output={idx}.{job_prefix_arg}.out\n"
            f"#SBATCH --error={idx}.{job_prefix_arg}.err\n"
            f"\nsource activate fastqc-env\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = fastqc_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


def create_unicycler_jobs(
    num_threads_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    isolate_list_arg,
    chunk_size_arg,
    fastq_ending_arg,
    input_dir_arg,
    output_dir_arg,
    job_prefix_arg,
    jobs_dir_arg,
    mode_arg,
    depth_filter_arg,
    min_contig_len_arg,
):

    """Function to create unicycler assembly jobs. Uses unicycler v0.5.0."""

    print(f'\nnum_threads_arg = "{num_threads_arg}"')
    print(f'memmax_arg = "{memmax_arg}"')
    print(f'maxtime_arg = "{maxtime_arg}"')
    print(f'partitions_arg = "{partitions_arg}"')
    print(f'isolate_list_arg = "{isolate_list_arg}"')
    print(f'chunk_size_arg = "{chunk_size_arg}"')
    print(f'fastq_ending_arg = "{fastq_ending_arg}"')
    print(f'input_dir_arg = "{input_dir_arg}"')
    print(f'output_dir_arg = "{output_dir_arg}"')
    print(f'job_prefix_arg = "{job_prefix_arg}"')
    print(f'jobs_dir_arg = "{jobs_dir_arg}"')
    print(f'mod_arg = "{mode_arg}"')
    print(f'depth_filter_arg = "{depth_filter_arg}"')
    print(f'min_contig_len_arg = "{min_contig_len_arg}"')

    # Read in isolate list:
    with open(isolate_list_arg, "r") as isolate_list_file:
        isolate_list = [isolate.strip() for isolate in isolate_list_file]

    # Create fastqc commands:
    unicycler_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(input_dir_arg, f"{isolate}_1{fastq_ending_arg}")
        input_r2 = os.path.join(input_dir_arg, f"{isolate}_2{fastq_ending_arg}")
        output_dir = os.path.join(output_dir_arg, isolate)

        unicycler_cmd = (
            f"unicycler -1 {input_r1} -2 {input_r2} -o {output_dir} "
            f"-t {num_threads_arg} --depth_filter {depth_filter_arg} "
            f"--min_fasta_len {min_contig_len_arg} --mode {mode_arg}\n"
        )

        unicycler_commands[isolate] = unicycler_cmd

    # Create chunked job files:
    chunks = [
        isolate_list[i : i + chunk_size_arg]
        for i in range(0, len(isolate_list), chunk_size_arg)
    ]

    for idx, chunk in enumerate(chunks):
        job_script = os.path.join(jobs_dir_arg, f"{idx}.{job_prefix_arg}.slurm")

        # Create job header:
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={maxtime_arg}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={num_threads_arg}\n"
            f"#SBATCH --partition={partitions_arg}\n"
            f"#SBATCH --mem={memmax_arg}\n"
            f"#SBATCH --output={idx}.{job_prefix_arg}.out\n"
            f"#SBATCH --error={idx}.{job_prefix_arg}.err\n"
            f"\nsource activate unicycler-env\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = unicycler_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


# --------------------------------------------------
if __name__ == "__main__":
    main()
