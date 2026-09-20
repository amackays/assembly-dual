Ava Mackay-Smith, genome assembly pipeline 2025

NOTE: this pipeline is set up for PARTIALLY PHASED assemblies only. If you are looking for a reference/alternate assembly pipeline (i.e. hifiasm --primary), see the [other repo](https://github.com/amackays/assembly-pri), which will purge and validate only the primary assembly of a Hifiasm assembly run.

# Setting up this Pipeline

First, after cloning this repo, you will need a conda environment that recognizes Snakemake (and, optionally, SLURM for job submission on the computing cluster).
You can use the .yaml file provided in this repo to make that environment: `conda env create -f SETUP.yaml`

Second, modify the permissions of the `SETUP-pull.sh` file with something like `chmod +x SETUP-pull.sh` and then run that bash script with `./SETUP-pull.sh` 
to download the larger attachments and Singularity images required to run the pipeline.

If your PacBio reads do not have informative file names, rename them to be sample names or similar, e.g. `cp 123445.bam.gz rosina.bam.gz` and store a SAFE COPY OF THIS FILE SOMEWHERE ELSE. NOW.

# Setup Structure of this Repo

This directory has the following subdirectories:
  - the directory `envs` contains specs for all conda environments required by rules inside the Snakefile.
  - the directory `scripts` contains all bash/python/R scripts called in the Snakefile.
  - the directory `slurm_general` houses a config file to make the Snakefile use SLURM-specific memory allocations.
  - optionally, the directory containing .bam Revio reads, e.g. `hifi_reads/`, can be located in this directory or somewhere else (the path must be included in the `config-ref.yaml` file, or whatever config file is sourced).

Runtime-generated directories:
  - the `logs` directory contains rule-specific logs which are created at runtime.
  - `analysis` may not exist yet, but is where all draft genomes and stats will be written.
  - a directory named `busco_downloads` will similarly appear after the pipeline is run.

The following files are also required to be present at runtime:
  - `config_ref.yaml` provides user-adjustable parameters for running the workflow, e.g. the path to .bam read files. The name of this file can be changed, but would also need to be changed in the Snakefile header.
  - the `Snakefile` contains all the rule chains for running the workflow.

Optional files and folders:
  - at runtime, if the SLURM profile is used, there will be a top-level log of rule order and outcome written to a file with syntax like `assembly-test.jobnumber.out`.
  - I keep a file at the top of the directory with the original .bam file names and the samples they correspond to.
  - the script `ref-slurm-wrap.sh` can be modified for use as a job submission wrapper for SLURM.

# To Run the Pipeline

You will first need to modify the config file you are using (e.g. `config-ref.yaml`) to contain paths to your reads, etc. The script `ref-slurm-wrap.sh` has a basic usage example of how a job can be run from the top level of this directory.

If launching interactively, I still strongly recommend launching Snakemake with SLURM for memory allocation. `snakemake --workflow-profile slurm_general`

If you do not use SLURM, make sure to request conda env usage: `snakemake --use-conda` and add whatever other parameters you would like.

# Troubleshooting

If you are using a computing cluster, make sure that your python install location is in your active conda env or the env of the rule in question, rather than installed through your cluster's Anaconda install.
This fixed many of my problems.

Other known bugs: if a job fails or the snakemake run is cut short (a user aborts a run, etc). then sometimes the conda environments built from the manifests in envs/ will fail with an error that Snakemake can't find an installed tool in that env. This can be solved by removing the built env and reinstalling it/forcing a reinstall.

Some of the late-stage tools are buggy (e.g. MitoFinder and NCBI-FCS). MitoFinder should be configured that it continues the run even if no mitochondrial contig is identified, but I haven't extensively validated this, partially because the MitoHiFi tool would be better for this job - I've left MitoFinder in this script to make it follow the workflow of the other [reference genome assembly pipeline](https://github.com/amackays/assembly-pri). If NCBI-FCS needs to be skipped, there's a separate Snakemake file in this directory called Snakefile_nofcs.

Feel free to modify the slurm_general config to work for your genomes (e.g. runtime, restart-times, etc).

[Notes on GitHub markdown formatting for READMEs]: (https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
