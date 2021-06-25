import os.path
import sys
from asm_hifi import loadPreCmd
from asm_hifi.snakemake_helper import DEFAULT_HPC_CONFIG_FILE
from asm_hifi import HpcConfig

def get_sample(wc):
	return wc.sample

def get_trimmed_reads(wc):
	return reads[wc.sample]

def get_read_file(wc):
	return reads_lookup[wc.sample][wc.read_basename]

def get_assembly(wc):
	return assembly_lookup[wc.sample][wc.assembly_basename]

def get_assem_bname(wc):
	return wc.assembly_basename

def get_alignments(wc):
	return alignments[wc.sample][wc.assembly_basename]

HPC_CONFIG = HpcConfig(DEFAULT_HPC_CONFIG_FILE)
TARGETS = []
assemblies = {}
reads = {}
reads_lookup = {}
assembly_lookup = {}
infile = open(config['sample_sheet'])
for line in infile:
	if line.startswith("#"):
		continue
	fields = line.rstrip().split(",")
	num_fields = len(fields)
	if num_fields != 3:
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample sheet should contain exactly 3 comma delimted fields: sample_name,reads_fofn,assemblies_fofn")
	
	sample_name = fields[0]
	reads_fof_name = fields[1]
	assem_fof_name = fields[2]

	if not os.path.isfile(reads_fof_name):
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": reads fofn " + reads_fof_name + "does not exist")

	if not os.path.isfile(assem_fof_name):
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": assembly fofn " + assem_fof_name + "does not exist")

	if sample_name in reads.keys():
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample names must be unique")

	reads[sample_name] = []
	reads_lookup[sample_name] = {}
	assembly_lookup[sample_name] = {}
	assemblies[sample_name] = []

	reads_fofn = open(reads_fof_name)
	for l in reads_fofn:
		reads_file = l.rstrip()
		if not os.path.isfile(reads_file):	
			sys.exit("ERROR in reads fofn" + reads_fof_name + " read file " + reads_file + " does not exist")
		reads[sample_name].append(reads_file)
		rbname = os.path.basename(reads_file)
		if rbname in reads_lookup[sample_name].keys():
			sys.exit("ERROR in reads fofn" + reads_fof_name + " read file basenames must be unique")
		reads_lookup[sample_name][rbname] = reads_file
	
	assem_fofn = open(assem_fof_name)
	for assem_file in assem_fofn:
		assem_file = assem_file.rstrip()
		if not os.path.isfile(assem_file):
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " assembly file " + assem_file + " does not exist")
		assemblies[sample_name].append(assem_file.rstrip())
		abname = os.path.basename(assem_file)
		if abname in assembly_lookup[sample_name].keys():
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " assembly file basenames must be unique")
		assembly_lookup[sample_name][abname] = assem_file

for sample_name in reads.keys():
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".ccs.yak"))

alignments = {}
for sample_name in assemblies.keys():
	alignments[sample_name] = {}
	for assembly in assemblies[sample_name]:
		assem_bname = os.path.basename(assembly)
		alignments[sample_name][assem_bname] = []
		if "all.contigs" in assem_bname:
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"yak.qv.out"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.gaps.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.ev.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pchlst.ctg.bed"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"kat22.comp.stats"))
		for reads_file in reads[sample_name]:
			read_bname = os.path.basename(reads_file)
			if "all.contigs" in assem_bname: 
				TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))
				alignments[sample_name][assem_bname].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))
				
localrules: all, asset_step1, asset_step2, asset_step3, asset_step4

rule all:
	input: TARGETS

rule yak_hash:
	input:
		get_trimmed_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}" + ".ccs.yak")
	params:
		load = loadPreCmd("source yak-0.1_CBG")
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("yak_hash") * attempt
	shell:
		"{params.load} yak count -t {threads} -o {output} {input} || echo 'YAK FAILED' > {output}"
		

rule yak_qv:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}" + ".ccs.yak"),
		get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","yak.qv.out")
	params:
		load = loadPreCmd("source yak-0.1_CBG")
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("yak_qv") * attempt
	shell:
		"{params.load} yak qv -t {threads} {input[0]} {input[1]} > {output} || echo 'YAK FAILED' > {output}"

rule minimap2paf:
	input:
		get_read_file,
		get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","{read_basename}" + ".paf")
	params:
		load = loadPreCmd("source minimap2-2.11")
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minimap2paf") * attempt
	shell:
		"{params.load} minimap2 -t {threads} -x asm20 {input[1]} {input[0]} > {output} || echo 'MINIMAP FAILED' > {output}"

rule asset_step1:
	input:
		assembly = get_assembly,
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.gaps.bed")
	params:
		load = loadPreCmd("source asset-1.0.0_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		1
	shell:
		"{params.load} cd {params.outdir}/{params.sample}/QC/{params.assem_bname} && detgaps {input.assembly} > {output} || echo 'ASSET FAILED' > {output}"

rule asset_step2:
	input:
		alignments = get_alignments,
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.pb.bed")
	params:
		load = loadPreCmd("source asset-1.0.0_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		1
	shell:
		"{params.load} cd {params.outdir}/{params.sample}/QC/{params.assem_bname} && ast_pb {input.alignments} > {output} || echo 'ASSET FAILED' > {output}"

rule asset_step3:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.gaps.bed"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.pb.bed")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.ev.pb.bed")
	params:
		load = loadPreCmd("source asset-1.0.0_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		1
	shell:
		"{params.load} cd {params.outdir}/{params.sample}/QC/{params.assem_bname} && acc {input[0]} {input[1]} > {output} || echo 'ASSET FAILED' > {output}"

rule asset_step4:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.gaps.bed"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.ev.pb.bed")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.pchlst.ctg.bed")
	params:
		load = loadPreCmd("source asset-1.0.0_CBG"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		1
	shell:
		"{params.load} cd {params.outdir}/{params.sample}/QC/{params.assem_bname} && pchlst -c {input[0]} {input[1]} > {output} || echo 'ASSET FAILED' > {output}"

rule kat:
	input:
		reads = get_trimmed_reads,
		assembly = get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","kat22.comp.stats")
	params:
		load = loadPreCmd("source kat-dev"),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("kat") * attempt
	shell:
		"{params.load} kat comp -t {threads} -m 22 -o {params.outdir}/{params.sample}/QC/{params.assem_bname}/kat22.comp '{input.reads}' {input.assembly} || echo 'KAT FAILED' > {output}"
