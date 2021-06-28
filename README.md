# DToL_asm_hifi
Assembly (Hifi) pipeline

## V1.1 includes 2 modules.

An Assembly Module does trimming of smrtbell adapters with cutadapt and assembly using hifiasm.
A QC module produces various basic assembly QC metrics on the assemblies produced.

## V1.2 includes the same 2 modules.

And adds a summary report summarising assembly metrics and QC metrics across assemblies.

## The branch 'feature_10x_functionality'
(Which is yet to become a release) currently adds one module which adds the functionality to polish assemblies with 10X data.
It is intended to also add in this branch a further module for 10X scaffolding / assembly breaking.
The polishing module is working code which has been tested on both Marigold and Tilapia genomes.

# Installation

## Python Dependencies

Requires python >=3.4

These dependencies should be installed automatically if not present.

Snakemake >= 4.4.0 
drmaa (for hpc environments)
matplotlib

## Installation from github

git clone https://github.com/EI-CoreBioinformatics/DToL_asm_hifi

cd DToL_asm_hifi

python setup.py bdist_wheel

pip install -I --prefix=/required_build_location dist/*whl

## On the EI cluster

source git-2.7.1

source python_miniconda-4.5.4_py3.6_gl

git clone https://github.com/EI-CoreBioinformatics/DToL_asm_hifi

cd DToL_asm_hifi/asm_hifi

python setup.py bdist_wheel

pip install -I --prefix=/required_build_location dist/*whl

## Put asm_hifi in your path on the EI cluster

export PYTHONPATH=/required_build_location/lib/python3.6/site-packages

echo PYTHONPATH=$PYTHONPATH

export PATH=/required_build_location/bin/:$PATH

# External Dependencies.

YAK

https://github.com/lh3/yak

hifiasm

https://github.com/chhylp123/hifiasm

cutadapt

https://cutadapt.readthedocs.io/en/stable/

minimap2

https://lh3.github.io/minimap2/minimap2.html

asset

https://github.com/dfguan/asset

KAT

https://kat.readthedocs.io/en/latest/walkthrough.html

bbmap

https://sourceforge.net/projects/bbmap/

# Config

The yaml configuration file asm_hifi/asm_hifi/etc/asm_hifi.config.yaml 
Currently includes parameter settings for cutadapt for trimming SMRT seqeunce adapters. It is not recommended to alter these settings.

The config file also includes instructions for how to source the various dependencies on your HPC.
If you are running on the EI cluster these can (and should) be left alone.

load:
        cutadapt: "source cutadapt-3.2_CBG"
        hifiasm: "source hifiasm-0.12"
        yak: "source yak-0.1_CBG"
        minimap2: "source minimap2-2.11"
        asset: "source asset-1.0.0_CBG"
        kat: "source kat-dev"
        bbmap: "source bbmap-38.06"
        
Otherwise replace these with the command to put thse dependencies in your path on your own system.
In a workstation / server enviroment where these are in your path just pass the empty string ""

The hpc config asm_hifi/asm_hifi/etc/hpc_config.json defines the resources you wish to allocate to the various steps in the pipeline.
This config SHOULD be edited to allocate sensible resources for the genome(s) you wish to assemble.
The defaults should be sufficient for a fairly small genome (they are tested with e.coli and should be sufficient for arabidopsis).
The Marigold genome is a 2.8gb tetraploid and assembles in hifiasm with 128gb of RAM. In under 48 hours using 32 cores.
Please see the hifiasm github for benchmarking on genomes of various sizes.
https://github.com/chhylp123/hifiasm

# Usage

## Module 1:

asm_hifi Assemble --outdir /output/directory sample_sheet

example sample sheet:

ecoli,4.6M,2,SRR10971019_subreads.fastq.gz

Running Module 1 creates the sample sheet for module 2.

## Module 2:

asm_hifi QC --outdir /output/directory sample_sheet

example sample sheet:

ecoli,/Module1_outdir/ecoli/ecoli.trimmed.reads.list,Module1_outdir/ecoli/ecoli.assem.outputs.list

## Module 3 (If 10x branch is checked out):

asm_hifi Polish -c 200 --outdir /output/directory sample_sheet

-c paramter controls max number of cluster cores to use (here this is effectively limiting the number of jobs)
It is important to set this. If not the pipeline will fill the cluster with freebayes jobs on an assembly with many contigs. 

Example sample sheet:

marigold	Assembly_all.contigs.fasta	/Reads/10X/R0073-S0006_Marigold_A10232_1_1_H32LTDRXY,/Reads/10X/R0073-S0006_Marigold_A10232_2_1_H32LTDRXY,/Reads/10X/R0073-S0006_Marigold_A10232_3_1_H32LTDRXY4

