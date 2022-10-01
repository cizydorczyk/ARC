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
        description="Create FastQC, Trimmomatic, Unicycler, SKESA, or Bakta jobs for ARC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # ----------------------------------------------
    # General options
    # ----------------------------------------------
    parser.add_argument(
        "-an",
        "--analysis",
        help="Analysis to run, one of: fastqc, trimmomatic, unicycler, skesa, or bakta. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-nt",
        "--num_threads",
        help="Number of threads to use for job. Set based on partitions desired. Default=8 appropriate for partitions single,lattice.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-mx",
        "--memmax",
        help="Max memory to use for job in MB; default=all available (0). Required.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-mt",
        "--maxtime",
        help="max time for job to run. Required.",
        type=str,
    )

    parser.add_argument(
        "-p",
        "--partitions",
        help='ARC partitions to use, comma-separated list; default="single,lattice". Required.',
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-i",
        "--isolate_list",
        help="path to isolate list on local computer; one isolate per line in list. Required.",
        metavar="file",
    )

    parser.add_argument(
        "-cs",
        "--chunk_size",
        help="number of isolates (chunk) to analyze in one job file. Required. Recommended size for Unicycler = 1.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-fe",
        "--fastq_ending",
        help="fastq file ending, one of: .fastq, .fq, .fastq.gz, or .fq.gz. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-pe",
        "--pair_ending",
        help="how forward/reverse read pairs are specified. One of _1,_2 or _R1,_R2.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-d",
        "--input_dir",
        help="input directory on ARC containing fastq files to be trimmed. Default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="output directory on ARC where trimmed files will be placed. Default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-j",
        "--jobs_dir",
        help="directory (local) to write job files to. default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-pre",
        "--prefix",
        help="job name prefix for job scripts, .out, and .err files. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-e",
        "--env",
        help="Conda environment on ARC that has desired analysis software installed. Required.",
        type=str,
        metavar="str",
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
        help="Minimum contig length Unicycler should produce. Required for Unicycler & Bakta.",
        type=int,
        metavar="int",
    )

    # ----------------------------------------------
    # Bakta options
    # ----------------------------------------------

    parser.add_argument(
        "-db",
        "--bakta_db",
        help="Path to Bakta db. Must already be downloaded & AMRFinderPlus portion initialized. Required for Bakta.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-g",
        "--genus",
        help="Genus for annotation. Required for Bakta.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-s",
        "--species",
        help="Species for annotation. Required for Bakta.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-as",
        "--assembly_dir",
        help="directory with assemblies for annotation. Required for Bakta.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "--assembly_suffix",
        help="Suffix of assembly filenames. E.g. if assembly file is JMB00961a-skesa-assembly.fasta, the suffix is 'skesa-assembly.fasta'. Leave out the first '-'! It will be placed automatically by the script. Must be a dash (for now)!",
        type=str,
        metavar="str",
    )

    return parser.parse_args()


# --------------------------------------------------


def main():
    """Create a variety of ARC jobs."""

    args = get_args()

    checkAnalysis(args.analysis)
    checkAnalysisOptions(args.analysis, args)

    if args.analysis == "trimmomatic":

        create_trimmomatic_jobs(
            num_threads_arg=args.num_threads,
            memmax_arg=args.memmax,
            maxtime_arg=args.maxtime,
            partitions_arg=args.partitions,
            isolate_list_arg=args.isolate_list,
            chunk_size_arg=args.chunk_size,
            fastq_ending_arg=args.fastq_ending,
            input_dir_arg=args.input_dir,
            output_dir_arg=args.output_dir,
            job_prefix_arg=args.job_prefix,
            adapters_file_arg=args.adapters_file,
            min_len_arg=args.min_len,
            min_qual_arg=args.min_qual,
            read_len_arg=args.read_len,
            jobs_dir_arg=args.jobs_dir,
            pair_ending_arg=args.pair_ending.split(","),
            env_arg=args.env,
        )

    if args.analysis == "fastqc":

        create_fastqc_jobs(
            num_threads_arg=args.num_threads,
            memmax_arg=args.memmax,
            maxtime_arg=args.maxtime,
            partitions_arg=args.partitions,
            isolate_list_arg=args.isolate_list,
            chunk_size_arg=args.chunk_size,
            fastq_ending_arg=args.fastq_ending,
            input_dir_arg=args.input_dir,
            output_dir_arg=args.output_dir,
            job_prefix_arg=args.prefix,
            jobs_dir_arg=args.jobs_dir,
            pair_ending_arg=args.pair_ending.split(","),
            env_arg=args.env,
        )

    if args.analysis == "unicycler":

        create_unicycler_jobs(
            num_threads_arg=args.num_threads,
            memmax_arg=args.memmax,
            maxtime_arg=args.maxtime,
            partitions_arg=args.partitions,
            isolate_list_arg=args.isolate_list,
            chunk_size_arg=args.chunk_size,
            fastq_ending_arg=args.fastq_ending,
            input_dir_arg=args.input_dir,
            output_dir_arg=args.output_dir,
            job_prefix_arg=args.prefix,
            jobs_dir_arg=args.jobs_dir,
            mode_arg=args.mode,
            depth_filter_arg=args.depth_filter,
            min_contig_len_arg=args.min_contig_len,
            pair_ending_arg=args.pair_ending.split(","),
            env_arg=args.env,
        )

    if args.analysis == "skesa":

        create_skesa_jobs(
            num_threads_arg=args.num_threads,
            memmax_arg=args.memmax,
            maxtime_arg=args.maxtime,
            partitions_arg=args.partitions,
            isolate_list_arg=args.isolate_list,
            chunk_size_arg=args.chunk_size,
            fastq_ending_arg=args.fastq_ending,
            input_dir_arg=args.input_dir,
            output_dir_arg=args.output_dir,
            job_prefix_arg=args.prefix,
            jobs_dir_arg=args.jobs_dir,
            pair_ending_arg=args.pair_ending.split(","),
            env_arg=args.env,
        )

    if args.analysis == "bakta":

        create_bakta_jobs(
            num_threads_arg=args.num_threads,
            memmax_arg=args.memmax,
            maxtime_arg=args.maxtime,
            partitions_arg=args.partitions,
            isolate_list_arg=args.isolate_list,
            chunk_size_arg=args.chunk_size,
            output_dir_arg=args.output_dir,
            job_prefix_arg=args.prefix,
            jobs_dir_arg=args.jobs_dir,
            env_arg=args.env,
            baktadb_arg=args.bakta_db,
            genus_arg=args.genus,
            species_arg=args.species,
            assembly_arg=args.assembly_dir,
            min_contig_len_arg=args.min_contig_len,
            assembly_suffix=args.assembly_suffix,
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
    pair_ending_arg,
    env_arg,
):
    """Function to create Trimmomatic jobs"""

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
    print(f'env_arg = "{env_arg}"')

    # Read in isolate list:
    isolate_list = readIsolateList(isolate_list_arg)

    # Create trimmomatic commands:
    trimmomatic_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        input_r2 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )
        output_r1 = os.path.join(
            output_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        output_r2 = os.path.join(
            output_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )
        output_r1_unpaired = os.path.join(
            output_dir_arg, f"{isolate}_u{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        output_r2_unpaired = os.path.join(
            output_dir_arg, f"{isolate}_u{pair_ending_arg[1]}{fastq_ending_arg}"
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
            f"\nsource activate {env_arg}\n\n"
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
    pair_ending_arg,
    env_arg,
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
    print(f'env_arg = "{env_arg}"')

    # Read in isolate list:
    isolate_list = readIsolateList(isolate_list_arg)

    # Create fastqc commands:
    fastqc_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        input_r2 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )

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
            f"\nsource activate {env_arg}\n\n"
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
    pair_ending_arg,
    env_arg,
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
    print(f'env_arg = "{env_arg}"')

    # Read in isolate list:
    isolate_list = readIsolateList(isolate_list_arg)

    # Create fastqc commands:
    unicycler_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        input_r2 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )
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
            f"\nsource activate {env_arg}\n\n"
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


def create_skesa_jobs(
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
    pair_ending_arg,
    env_arg,
):

    """Function to create SKESA assembly jobs. Uses SKESA v2.4.0."""

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
    print(f'env_arg = "{env_arg}"')

    # Read in isolate list:
    isolate_list = readIsolateList(isolate_list_arg)

    # Create fastqc commands:
    skesa_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        input_r2 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )
        output_assembly = os.path.join(
            output_dir_arg, f"{isolate}-skesa-assembly.fasta"
        )

        skesa_cmd = f"skesa --reads {input_r1},{input_r2} --cores {num_threads_arg} --contigs_out {output_assembly}"

        skesa_commands[isolate] = skesa_cmd

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
            f"\nsource activate {env_arg}\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = skesa_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


def create_bakta_jobs(
    num_threads_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    isolate_list_arg,
    chunk_size_arg,
    output_dir_arg,
    job_prefix_arg,
    jobs_dir_arg,
    env_arg,
    baktadb_arg,
    genus_arg,
    species_arg,
    assembly_arg,
    min_contig_len_arg,
    assembly_suffix,
):
    """Function to create Bakta jobs for ARC."""

    # Read in isolate list:
    isolate_list = readIsolateList(isolate_list_arg)

    # Create bakta commands:
    bakta_commands = {}
    for isolate in isolate_list:
        assembly = os.path.join(assembly_arg, f"{isolate}-{assembly_suffix}")
        output_dir = os.path.join(output_dir_arg, isolate)

        bakta_cmd = (
            f"bakta --db {baktadb_arg} --min-contig-len {min_contig_len_arg} --prefix {isolate} --output {output_dir} "
            f"--genus {genus_arg} --species {species_arg} --strain {isolate} --translation-table 11 --compliant "
            f"--threads {num_threads_arg} {assembly}"
        )

        bakta_commands[isolate] = bakta_cmd

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
            f"\nsource activate {env_arg}\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = bakta_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


def checkOption(option_name, option_arg=None):
    if option_arg is None:
        print(
            f"\nError: option {option_name} is required for your selected analysis type. "
            f"Please provide {option_name} and run again. Exiting.\n"
        )
        sys.exit()


def checkAnalysis(analysis_type):
    if analysis_type is None or analysis_type not in [
        "fastqc",
        "trimmomatic",
        "unicycler",
        "skesa",
        "bakta",
    ]:
        print(
            f"\nError: --analysis not specified or set to invalid type. Please select from fastqc, trimmomatic, "
            f"unicycler, or skesa. Exiting.\n"
        )
        sys.exit()


def readIsolateList(isolate_list_file):
    # Read in isolate list:
    with open(isolate_list_file, "r") as isolate_list_file:
        isolate_list = [isolate.strip() for isolate in isolate_list_file]

    return isolate_list


def checkAnalysisOptions(analysis_type, parsed_args):
    common_options_list = [
        ("--num_threads", parsed_args.num_threads),
        ("--memmax", parsed_args.memmax),
        ("--maxtime", parsed_args.maxtime),
        ("--partitions", parsed_args.partitions),
        ("--isolate_list", parsed_args.isolate_list),
        ("--chunk_size", parsed_args.chunk_size),
        ("--fastq_ending", parsed_args.fastq_ending),
        ("--pair_ending", parsed_args.pair_ending),
        ("--input_dir", parsed_args.input_dir),
        ("--output_dir", parsed_args.output_dir),
        ("--jobs_dir", parsed_args.jobs_dir),
        ("--prefix", parsed_args.prefix),
        ("--env", parsed_args.env),
    ]

    analysis_specific_options_dict = {
        "trimmomatic": common_options_list
        + [
            ("--adapters_file", parsed_args.adapters_file),
            ("--min_len", parsed_args.min_len),
            ("--read_len", parsed_args.read_len),
            ("--min_qual", parsed_args.min_qual),
        ],
        "fastqc": common_options_list,
        "unicycler": common_options_list
        + [
            ("--mode", parsed_args.mode),
            ("--depth_filter", parsed_args.depth_filter),
            ("--min_contig_len", parsed_args.min_contig_len),
        ],
        "skesa": common_options_list,
        "bakta": [
            ("--num_threads", parsed_args.num_threads),
            ("--memmax", parsed_args.memmax),
            ("--maxtime", parsed_args.maxtime),
            ("--partitions", parsed_args.partitions),
            ("--isolate_list", parsed_args.isolate_list),
            ("--chunk_size", parsed_args.chunk_size),
            ("--output_dir", parsed_args.output_dir),
            ("--jobs_dir", parsed_args.jobs_dir),
            ("--prefix", parsed_args.prefix),
            ("--env", parsed_args.env),
            ("--genus", parsed_args.genus),
            ("--species", parsed_args.species),
            ("--bakta_db", parsed_args.bakta_db),
            ("--assembly_dir", parsed_args.assembly_dir),
            ("--assembly_suffix", parsed_args.assembly_suffix),
            ("--min_contig_len", parsed_args.min_contig_len),
        ],
    }

    for option in analysis_specific_options_dict[analysis_type]:
        checkOption(option[0], option[1])


# --------------------------------------------------
if __name__ == "__main__":
    main()
