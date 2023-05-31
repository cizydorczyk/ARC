# Scripts for creating jobs on ARC.

This repository contains scripts (only one currently) for creating **FastQC**, **Trimmomatic**, **Unicycler** (assembly), **SKESA** (assembly), and **Bakta** (annotation) jobs for ARC. These can then be uploaded to ARC and submitted to run the desired analyses.

See the [Wiki](https://github.com/cizydorczyk/ARC/wiki) for detailed guides.

## Quick Guide
This script will make job scripts that can be uploaded to ARC and submitted. There is the option to specify multiple analyses per job script (the `--chunk_size` option), e.g. assemble two genomes per job script.  

Job scripts are named with an integer followed by a prefix and ending with the extension '.slurm'. For example, '0.trimmomatic.slurm'. The integer is meaningless; the order of samples in your isolate list (`--isolate_list` option) determines which sample is in which file. I.e., the first sample will be in 0.trimmomatic.slurm, the second in 1.trimmomatic.slurm, etc. (if chunk size is 1).  

### FastQ File Naming Conventions
FastQ files must include three parts: the sample name, the fastq pair ending, and the file extensions.  

For example, a file may be named 'Sample01_R1.fastq.gz'. In this case, the 'Sample01' is the sample name, the '\_R1' is the pair ending, and the '.fastq.gz' is the file extensions.

Sample names can include any common characters (including underscores and dashes but no spaces). Pair endings should be one of '\_R1' and '\_R2' or '\_1' and '\_2'. File extensions should be .fastq or .fq and can be optionally gzipped: .fastq.gz or .fq.gz.

The following are some examples of acceptable fastq file name formats:  
Sample01_R1.fastq  
Sample_01_1.fastq.gz  
Sample01_1.fq 
Sample_01_R1.fq.gz  

Note that the entirety of the sample name must be included in the isolate list (see options).  

## Example Commands


### Trimmomatic example command:

```
$ create_arc_fastqc_trimmomatic_unicycler_jobs.py --analysis trimmomatic --num_threads 12 --memmax 0 --maxtime 04:00:00 --partitions parallel --isolate_list pa-novogene-apr2023-isolate-list.txt --chunk_size 4 --fastq_ending .fastq.gz --pair_ending _1,_2 --input_dir /work/parkins_lab/project/conrad/all-pa/fastq-files/raw-fastq/pa-novogene-apr2023-raw-fastq/ --output_dir /work/parkins_lab/project/conrad/all-pa/fastq-files/trimmed-fastq/ --jobs_dir pa-novogene-apr2023-trim-jobs/ --prefix pangapr23 --env trimmomatic-env --adapters_file /home/conrad.izydorczyk/TruSeq3-PE-2.fa --min_len 31 --read_len 150 --min_qual 5
```

### Unicycler example command:

```
$ create_arc_fastqc_trimmomatic_unicycler_jobs.py --analysis unicycler --num_threads 8 --memmax 0 --maxtime 24:00:00 --partitions single,lattice --isolate_list ~/pa-novogene-jan2023/pa-novogene-jan2023-isolate-list.txt --chunk_size 1 --fastq_ending .fastq.gz --pair_ending _1,_2 --input_dir /work/parkins_lab/project/conrad/pa-novogene-jan2023/fastq-files/z6-normalized-fastq/ --output_dir /work/parkins_lab/project/conrad/pa-novogene-jan2023/de-novo-assemblies/z7-unicycler/ --jobs_dir z7-unicycler-jobs/ --prefix uni --env unicycler-env --mode normal --depth_filter 0.25 --min_contig_len 200
```

### SKESA example command:

```
$ create_arc_fastqc_trimmomatic_unicycler_jobs.py --analysis skesa --num_threads 8 --memmax 0 --maxtime 24:00:00 --partitions single,lattice --isolate_list ~/pa-shared-2017-2020/pa-les-isolate-list.txt --chunk_size 2 --fastq_ending .fastq.gz --pair_ending _1,_2 --input_dir /work/parkins_lab/project/conrad/pa-shared-2017-2020/fastq-files/o2-trimmed-fastq/ --output_dir /work/parkins_lab/project/conrad/pa-shared-2017-2020/de-novo-assemblies/o4-skesa/ --jobs_dir o4-skesa-jobs/ --prefix sksa --env skesa-env
```

### Bakta example command:

```
$ create_arc_fastqc_trimmomatic_unicycler_jobs.py --analysis bakta --num_threads 8 --memmax 0 --maxtime 04:00:00 --partitions single,lattice --isolate_list ~/pa-shared-2017-2020/pa-les-isolate-list.txt --chunk_size 1 --output_dir /work/parkins_lab/project/conrad/pa-shared-2017-2020/genome-assemblies/o5-bakta-unicycler/ --jobs_dir o5-bakta-unicycler-jobs/ --prefix unibkt --env bakta-env --genus Pseudomonas --species aeruginosa --bakta_db /home/conrad.izydorczyk/bakta-db-2022-08-29/db/ --assembly_dir /work/parkins_lab/project/conrad/pa-shared-2017-2020/de-novo-assemblies/o4-unicycler/ --assembly_suffix assembly.fasta --min_contig_len 200
```

## Option Explanations
There are many options, some universal to all analyses and some specific to a given tool.

### General (universal) options:
`--analysis`: Which analysis to run. One of 'fastqc', 'trimmomatic', 'unicycler', 'skesa', or 'bakta'.  
`--num_threads`: Number of threads to use for both the job request on ARC and for the tool being run. This should correspond to # of threads available on desired partitions on ARC.  
`--memmax`: Max amount of memory for job. Set max based on ARC partition, or set to '0' to use all available memory in requested partition.  
`--maxtime`: Max time for job to run. Set based on max times available on ARC partitions & how long you expect an analysis to take. FastQC/Trimmomatic are quick, so if you are running a job script *per sample*, setting this low is okay (15min-1hr). If you are assembling genomes or annotating, suggestions are 2-4hr for Unicycer, 1-2hr for SKESA, and 4-5hr for Bakta. These times ensure that your job will finish; these tools generally run *much faster* than this.  
`--partitions`: Partitions you wish to use on ARC. Comma separated list like so: `--partitions single,lattice`.  
`--isolate_list`: Plain text file with one sample per line.  
`--chunk_size`: Number of samples to process *per job*. For FastQC/Trimmoamtic, can set to 5-15. For Unicycler/Bakta, set to 1. For SKESA, can set to 1-4. Adjust `--maxtime` accordingly.  
`--input_dir`: Input directory with fastq files **on ARC**. Must follow naming conventions (see below & `--fastq_ending` and `--pair_ending`). Used for all analyses except Bakta.  
`--output_dir`: Output directory for analysis **on ARC**. Individual sample directories will be created as required.  
`--fastq_ending`: File extensions after \_R1 or \_1 in your fastq names. E.g. Sample-1_R1.fastq.gz = '.fastq.gz' (note the '.'s!). Used for all analyses except Bakta.  
`--pair_ending`: How your forward/reverse fastq files are differentiated. Generally, one of '\_R1,\_R2' or '\_1,\_2'. Comma separated, no spaces. Used for all analyses except Bakta.  
`--jobs_dir`: Directory for jobs scripts **on local computer**. Upload this directory and submit job scripts on ARC.  
`--prefix`: Prefix to append to job script names. E.g., for Fastqc, could be 'fqc'.  
`--env`: Name of conda environment in which tools required for desired analysis are installed. E.g., for an environment with Bakta version 1.7 installed in a conda environment named 'bakta-env', this would be `--env bakta-env`.  















