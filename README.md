# DToL_asm_hifi
Assembly (Hifi) pipeline

## V1.1 includes 2 modules.

An Assembly Module does trimming of smrtbell adapters with cutadapt and assembly using hifiasm.
A QC module produces various basic assembly QC metrics on the assemblies produced.

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

source python_miniconda-4.5.4_py3.6_gl

export PYTHONPATH=/required_build_location/lib/python3.6/site-packages

echo PYTHONPATH=$PYTHONPATH

export PATH=/required_build_location/bin/:$PATH

# Usage

## Module 1:

asm_hifi Assemble --outdir /output/directory sample_sheet

## Module 2:

asm_hifi QC --outdir /output/directory sample_sheet

## Module 3 (If 10x branch is checked out):

asm_hifi Polish -c 200 --outdir /output/directory sample_sheet

-c paramter controls max number of cluster cores to use (here this is effectively limiting the number of jobs)
It is important to set this. If not the pipeline will fill the cluster with freebayes jobs on an assembly with many contigs. 

