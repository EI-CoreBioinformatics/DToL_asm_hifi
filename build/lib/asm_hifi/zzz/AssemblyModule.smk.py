import os.path
import sys
from asm_hifi import loadPreCmd
from asm_hifi.snakemake_helper import DEFAULT_HPC_CONFIG_FILE
from asm_hifi import HpcConfig

def get_sample(wc):
	return wc.sample

def get_purge_level(wc):
	return wc.purge_level

def get_read_file(wc):
	return reads[wc.sample][wc.read_bname]

def get_trimmed_reads(wc):
	return trimmed_reads[wc.sample]

def get_purge_arg(wc):
	if wc.purge_level == "L0":
		return "-l 0"
	elif wc.purge_level == "L1":
		return "-l 1"
	elif wc.purge_level == "L2":
		return "-l 2"

def get_sample_final_assemblies(wc):
	return assemblies_by_sample[wc.sample]

HPC_CONFIG = HpcConfig(DEFAULT_HPC_CONFIG_FILE)
TARGETS = []
assemblies = {}
reads = {}
trimmed_reads = {}
final_assemblies = []
assemblies_by_sample = {}
output_lists = []

#read sample sheet and populate snakemakes 'rule all'
infile = open(config['sample_sheet'])
for line in infile:
	if line.startswith("#"):
		continue
	fields = line.rstrip().split(",")
	num_fields = len(fields)
	if num_fields < 4:
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample sheet should contain at least 4 comma delimited fields\nsample_name,est_genome_size,ploidy,reads... [reads].. [reads]\ngenome_size and ploidy are optional and may be an empty field")
	sample_name = fields[0]
	#genome size is currently not used. TODO add code to check its of the form \dM||\dG
	genome_size = fields[1]
	if not fields[2].isdigit():
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": ploidy must be a positive integer")
	ploidy = int(fields[2])
	if sample_name in reads.keys():
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample names must be unique")
	reads[sample_name] = {}
	trimmed_reads[sample_name] = []
	assemblies_by_sample[sample_name] = []
	for i in range(3,num_fields):
		if not os.path.isfile(fields[i]):
			sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": reads file " + fields[i] + " does not exist.")
		bname = os.path.basename(fields[i])
		if bname in reads[sample_name].keys():
			sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": readfile basenames must be unique")
		reads[sample_name][bname] = fields[i]
		trimmed_reads[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"reads",bname + ".cutadapt.keep.fastq"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"reads",bname + ".cutadapt.keep.fastq"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"reads",bname + ".cutadapt.reject.fastq"))

	if ploidy == 1:
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.a_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.p_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.a_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.p_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.cleanup.done"))
		final_assemblies.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.p_ctg.fasta"))
		assemblies_by_sample[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L0.p_ctg.fasta"))

	else:
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.a_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.p_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.a_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.p_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.all.contigs.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.cleanup.done"))
		final_assemblies.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.p_ctg.fasta"))
		final_assemblies.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.all.contigs.fasta"))
		assemblies_by_sample[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.p_ctg.fasta"))
		assemblies_by_sample[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L1.all.contigs.fasta"))

		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.a_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.p_ctg.gfa"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.a_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.p_ctg.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.all.contigs.fasta"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.cleanup.done"))
		final_assemblies.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.p_ctg.fasta"))
		final_assemblies.append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.all.contigs.fasta"))
		assemblies_by_sample[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.p_ctg.fasta"))
		assemblies_by_sample[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"hifiasm",sample_name + ".L2.all.contigs.fasta"))

	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".assem.outputs.list"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".trimmed.reads.list"))
	output_lists.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".assem.outputs.list"))
	output_lists.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".trimmed.reads.list"))

TARGETS.append(os.path.join(config['output_base_dir'],"qc.sample.sheet"))

localrules: all, hifiasm_cleanup, create_fasta, list_assem, list_reads, list_outputs

rule all:
	input:
		TARGETS

rule cutadapt:
	input:
		get_read_file
	output:
		os.path.join(config['output_base_dir'],"{sample}","reads","{read_bname}" + ".cutadapt.keep.fastq"),
		os.path.join(config['output_base_dir'],"{sample}","reads","{read_bname}" + ".cutadapt.reject.fastq")
	params:
		load = loadPreCmd("source cutadapt-3.2_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		cutadapt_params = config['cutadapt_params']
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("cutadapt") * attempt
	shell:
		"{params.load} mkdir -p {params.outdir}/{params.sample}/reads && cutadapt --untrimmed-output {output[0]} {params.cutadapt_params} -o {output[1]} --rc {input}"

rule hifiasm:
	input:
		get_trimmed_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.a_ctg.gfa"),
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.p_ctg.gfa")
	params:
		outdir = config['output_base_dir'],
		sample = get_sample,
		purge_level = get_purge_level,
		purge_arg = get_purge_arg,
		load = loadPreCmd("source hifiasm-0.12")
	threads:
		8
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("hifiasm") * attempt
	shell:
		"{params.load} mkdir -p {params.outdir}/{params.sample}/hifiasm && hifiasm -t {threads} {params.purge_arg} -o {params.outdir}/{params.sample}/hifiasm/{params.sample}.{params.purge_level} {params.outdir}/{params.sample}/reads/*.cutadapt.keep.fastq"

if config['HIFIASM_KEEP_BIN_FILES']:

	rule hifiasm_cleanup:
		input:
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.a_ctg.gfa"),
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.p_ctg.gfa")
		output:
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.cleanup.done")
		params:
			outdir = config['output_base_dir'],
			sample = get_sample,
		threads:
			1
		shell:
			"touch {output}"
		
else:
	rule hifiasm_cleanup:
		input:
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.a_ctg.gfa"),
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.p_ctg.gfa")
		output:
			os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.cleanup.done")
		params:
			outdir = config['output_base_dir'],
			sample = get_sample,
			purge_level = get_purge_level,
		threads:
			1
		shell:
			"rm -f {params.outdir}/{params.sample}/hifiasm/{params.sample}.{params.purge_level}*.bin && touch {output}"

rule create_fasta:
	input:
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.cleanup.done"),
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.a_ctg.gfa"),
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.p_ctg.gfa")
	output:
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.a_ctg.fasta"),
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.p_ctg.fasta"),
		os.path.join(config['output_base_dir'],"{sample}","hifiasm","{sample}" + ".{purge_level}.all.contigs.fasta")
	params:
		outdir = config['output_base_dir'],
		sample = get_sample,
		awk_cmd = "awk '$1 ~/S/ {print \">\"$2\"\\n\"$3}'"
	shell:
		"cat {input[1]} | {params.awk_cmd} > {output[0]} && cat {input[2]} | {params.awk_cmd} > {output[1]} && cat {output[0]} {output[1]} > {output[2]}"

rule list_assem:
	input:
		get_sample_final_assemblies
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}" + ".assem.outputs.list")
	run:
		with open(output[0], "wt") as outfile:
			for input_path in input:
				outfile.write(input_path + "\n")
			
rule list_reads:
	input:
		get_trimmed_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}" + ".trimmed.reads.list")
	run:
		with open(output[0], "wt") as outfile:
			for input_path in input:
				outfile.write(input_path + "\n")

rule list_outputs:
	input:
		output_lists
	output:
		os.path.join(config['output_base_dir'],"qc.sample.sheet")
	run:
		filemap = {}
		for input_path in input:
			path_diff = os.path.relpath(input_path,config['output_base_dir'])
			sample = os.path.dirname(path_diff)
			if sample not in filemap.keys():
				filemap[sample] = {}
			if input_path.endswith(".outputs.list"):
				filemap[sample]['assem'] = input_path
			elif input_path.endswith(".reads.list"):
				filemap[sample]['reads'] = input_path

		with open(output[0], "wt") as outfile:
			for sample in filemap.keys():
				outfile.write(sample + "," + filemap[sample]['reads'] + "," + filemap[sample]['assem'] + "\n")
				
