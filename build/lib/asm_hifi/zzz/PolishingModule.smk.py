import os.path
import sys
from asm_hifi import loadPreCmd
from asm_hifi.snakemake_helper import DEFAULT_HPC_CONFIG_FILE
from asm_hifi import HpcConfig

def get_sample(wc):
	return wc.sample

def get_assembly_path(wc):
	return assemblies[wc.sample]

def get_assembly_basename(wc):
	return os.path.splitext(os.path.basename(assemblies[wc.sample]))[0]

def get_reads(wc):
	return reads[wc.sample]

def get_reads_string(wc):
	return os.path.dirname(reads[wc.sample].split(",")[0])

def get_lane_string(wc):
	lanelist = []
	for fastqdir in reads[wc.sample].split(","):
		lanelist.append(os.path.basename(fastqdir))
	return ",".join(lanelist)

localrules: all, mkdir_and_copy_asm, index_fasta, mkref, mark_contigs, list_bcf, vcf_index, bcftools_consensus
HPC_CONFIG = HpcConfig(DEFAULT_HPC_CONFIG_FILE)
TARGETS = []
assemblies = {}
reads = {}

infile = open(config['assembly_sheet'])
for line in infile:
	if line.startswith("#"):
		continue
	fields = line.rstrip().split("\t")
	sample_name = fields[0]
	assembly_path = fields[1]
	readsdir = fields[2]
	if "," in readsdir:
		read_dirs = readsdir.split(",")
		base_dir = os.path.dirname(read_dirs[0])
		subdirs = []
		for directory_name in read_dirs:
			if os.path.dirname(directory_name) != base_dir:
				sys.exit("ERROR IN SAMPLESHEET: All Fastq directrories should be sub directories of the same base directory")
			subdirs.append(os.path.basename(directory_name))
		reads[sample_name] = "--fastqs=" + base_dir + " --sample=" + ",".join(subdirs)
	else:
		reads[sample_name] = readsdir
	print(reads[sample_name])
	assemblies[sample_name] = assembly_path
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".fa"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".fa.fai"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".mkref.done"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name,"outs","possorted_bam.bam"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".bcf"))

	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".changes.vcf.gz"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".changes.vcf.gz.csi"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".polished.fasta"))

print(TARGETS)

rule all:
	input: TARGETS

rule mkdir_and_copy_asm:
	input:
		get_assembly_path
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa")
	params:
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	shell:
		"mkdir -p {params.outdir}/{params.sample} && cp {input} {output}"

rule index_fasta:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa")
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa.fai")
	params:
		load = loadPreCmd("source samtools-1.9_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	shell:
		"{params.load} samtools faidx {input}"

rule mkref:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa")
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.mkref.done")
	params:
		load = loadPreCmd("source longranger-2.1.2"),
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	shell:
		"{params.load} rm -rf {params.outdir}/{params.sample}/refdata-{params.sample} && mkdir -p {params.outdir}/{params.sample} && cd {params.outdir}/{params.sample} && longranger mkref {input} && touch {output}"

rule align:
	input:
		refdone = os.path.join(config['output_base_dir'],"{sample}","{sample}.mkref.done"),
		#reads = get_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}","outs","possorted_bam.bam"),
		os.path.join(config['output_base_dir'],"{sample}","{sample}","outs","summary.csv")
	params:
		load = loadPreCmd("source longranger-2.1.2"),
		assembly_basename = get_assembly_basename,
		outdir = config['output_base_dir'],
		sample = get_sample,
		read_string = get_reads,
		#lane_string = get_lane_string
	threads:
		32
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("align") * attempt
	shell:
		#"{params.load} cd {params.outdir}/{params.sample} && rm -rf {params.outdir}/{params.sample}/{params.sample} && longranger align --localcores={threads} --id={params.sample} --reference={params.outdir}/{params.sample}/refdata-{params.sample} --fastqs={params.read_string} --sample={params.lane_string}"
		"{params.load} cd {params.outdir}/{params.sample} && rm -rf {params.outdir}/{params.sample}/{params.sample} && longranger align --localcores={threads} --id={params.sample} --reference={params.outdir}/{params.sample}/refdata-{params.sample} {params.read_string}"

rule mark_contigs:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa.fai"),
		os.path.join(config['output_base_dir'],"{sample}","{sample}","outs","summary.csv"),
		os.path.join(config['output_base_dir'],"{sample}","{sample}","outs","possorted_bam.bam")
	output:
		dynamic(os.path.join(config['output_base_dir'],"{sample}","MarkContigs","contig_{i}.txt"))
	params:
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	run:
		mean_depth = 0
		linenum = 0
		with open(input[1],"r") as summary_csv:
			for line in summary_csv:
				linenum = linenum + 1
				fields = line.split(",")
				if linenum == 2:
					mean_depth = float(fields[16])
		max_depth = int(mean_depth * 5)
		with open(input[0],"r") as faindex:
			for line in faindex:
				fields = line.split("\t")
				contig_name = fields[0]
				contig_length = fields[1]
				with open(os.path.join(params.outdir,params.sample,"MarkContigs","contig_" + contig_name + ".txt"),"w") as outfile:
					outfile.write(contig_name + ":1-" + contig_length + "\t" + str(max_depth) + "\t" + contig_name + "\n") 
			
rule freebayes:
	input:
		os.path.join(config['output_base_dir'],"{sample}","MarkContigs","contig_{i}.txt"),
	output:
		os.path.join(config['output_base_dir'],"{sample}","BCF","contig_{i}.bcf")
	params:
		load = loadPreCmd("source freebayes-1.3.1","source bcftools-1.9"),
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	resources:
		 mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("freebayes") * attempt
	shell:
		"{params.load} paramarray=( $(cat {input[0]} ) ) && freebayes --bam {params.outdir}/{params.sample}/{params.sample}/outs/possorted_bam.bam --region=${{paramarray[0]}} --skip-coverage ${{paramarray[1]}} -f {params.outdir}/{params.sample}/{params.sample}.fa | bcftools view --no-version -Ou > {output}"

rule list_bcf:
	input:
		dynamic(os.path.join(config['output_base_dir'],"{sample}","BCF","contig_{i}.bcf"))
	output:
		os.path.join(config['output_base_dir'],"{sample}","BCF","concat.list")
	params:
		outdir = config['output_base_dir'],
		sample = get_sample
	shell:
		"echo {params.outdir}/{params.sample}/BCF/*.bcf | xargs ls > {output}"

rule bcf_concat:
	input:
		os.path.join(config['output_base_dir'],"{sample}","BCF","concat.list") 
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.bcf")
	params:
		load = loadPreCmd("source bcftools-1.9"),
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		8
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("bcf_concat") * attempt
	shell:
		"{params.load} bcftools concat -f {input} | bcftools view -Ou -e'type=\"ref\"' | bcftools norm -Ob -f {params.outdir}/{params.sample}/refdata-{params.sample}/fasta/genome.fa -o {output} --threads {threads}"

rule vcf_changes:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.bcf")
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.changes.vcf.gz")
	params:
		load = loadPreCmd("source bcftools-1.9")
	threads:
		8
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("vcf_changes") * attempt
	shell:
		"{params.load} bcftools view -i 'QUAL>1 && (GT=\"AA\" || GT=\"Aa\")' -Oz --threads={threads} {input} > {output}"

rule vcf_index:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.changes.vcf.gz")
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.changes.vcf.gz.csi")
	params:
		load = loadPreCmd("source bcftools-1.9")
	threads:
		1
	shell:
		"{params.load} bcftools index {input}"

rule bcftools_consensus:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.changes.vcf.gz"),
		os.path.join(config['output_base_dir'],"{sample}","{sample}.changes.vcf.gz.csi")
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.polished.fasta")
	params:
		load = loadPreCmd("source bcftools-1.9"),
		outdir = config['output_base_dir'],
		sample = get_sample
	threads:
		1
	shell:
		"{params.load} bcftools consensus -Hla -f {params.outdir}/{params.sample}/refdata-{params.sample}/fasta/genome.fa {input[0]} > {output}"
