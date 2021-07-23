import os.path
import sys
from asm_hifi import loadPreCmd
from asm_hifi.snakemake_helper import DEFAULT_HPC_CONFIG_FILE
from asm_hifi import HpcConfig
import pandas as pd

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

def get_busco_db(wc):
	return config["busco_dbs"][kingdoms[wc.sample]]

def get_results_files(wc):
	return results_files[wc.sample]

def extract_stats(statsfilename):
	with open(statsfilename,"r") as statsfile:
		for line in statsfile:
			if ":" in line:
				fields = line.rstrip().split("\t")
				if "Main genome contig total" in fields[0]:
					num_contigs = fields[1]
				elif "Main genome contig sequence total" in fields[0]:
					assem_size = fields[1].rstrip()
				elif "Main genome contig N/L50" in fields[0]:
					(N50,L50) = fields[1].split("/")
				elif "Main genome contig N/L90" in fields[0]:
					(N90,L90) = fields[1].split("/")
				elif "Number of scaffolds > 50 KB" in fields[0]:
					num_over50kb = fields[1]
				elif "% main genome in scaffolds > 50 KB" in fields[0]:
					perc_over50kb = fields[1]

	return ",".join([num_contigs,assem_size,L50,L90,num_over50kb,perc_over50kb])

def extract_qv(yakfilename):
	with open(yakfilename, "r") as yak:
		for line in yak:
			if line.startswith("QV"):
				fields = line.rstrip().split("\t")
				qv = fields[0]
			elif "YAK FAILED" in line:
				qv = "YAK FAILED"
	return qv

def extract_busco(buscofilename):
	i = 0
	with open(buscofilename, "r") as buscofile:
		for line in buscofile:
			i = i + 1
			if i == 8:
				busco = line.rstrip().replace("\t","").replace(","," ")
	return busco
			
				

HPC_CONFIG = HpcConfig(DEFAULT_HPC_CONFIG_FILE)
TARGETS = []
assemblies = {}
reads = {}
reads_lookup = {}
assembly_lookup = {}
#adapter_logs = []
kingdoms = {}
assembly_types = {}

infile = open(config['sample_sheet'])
for line in infile:
	if line.startswith("#"):
		continue
	fields = line.rstrip().split(",")
	num_fields = len(fields)
	if num_fields != 4:
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample sheet should contain exactly 4 comma delimted fields: sample_name,busco_clade,reads_fofn,assemblies_fofn")
	
	sample_name = fields[0]
	kingdom = fields[1]
	reads_fof_name = fields[2]
	assem_fof_name = fields[3]

	if not os.path.isfile(reads_fof_name):
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": reads fofn " + reads_fof_name + "does not exist")

	if not os.path.isfile(assem_fof_name):
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": assembly fofn " + assem_fof_name + "does not exist")

	if sample_name in reads.keys():
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample names must be unique")

	if sample_name in kingdoms.keys():
		sys.exit("ERROR in sample sheet" + config['sample_sheet'] + ": sample names must be unique")
	else:
		kingdoms[sample_name] = kingdom

	reads[sample_name] = []
	reads_lookup[sample_name] = {}
	assembly_lookup[sample_name] = {}
	assemblies[sample_name] = []
	if not sample_name in assembly_types.keys():
		assembly_types[sample_name] = {}

	reads_fofn = open(reads_fof_name)
	for l in reads_fofn:
		reads_file = l.rstrip()
		#adapter_logs.append(reads_file + ".log")
		if not os.path.isfile(reads_file):	
			sys.exit("ERROR in reads fofn" + reads_fof_name + " read file " + reads_file + " does not exist")
		reads[sample_name].append(reads_file)
		rbname = os.path.basename(reads_file)
		if rbname in reads_lookup[sample_name].keys():
			sys.exit("ERROR in reads fofn" + reads_fof_name + " read file basenames must be unique")
		reads_lookup[sample_name][rbname] = reads_file
	
	assem_fofn = open(assem_fof_name)
	for line in assem_fofn:
		fields = line.rstrip().split(",")
		if len(fields) != 2:
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " each entry must be comma delimited of the form assembly_type,assembly_file_path")
		assem_type = fields[0]
		assem_file = fields[1]
		if assem_type != "all" and assem_type != "primary":
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " assembly type must be 'primary' or 'all'")
		if not os.path.isfile(assem_file):
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " assembly file " + assem_file + " does not exist")
		assemblies[sample_name].append(assem_file.rstrip())
		abname = os.path.basename(assem_file)
		if abname in assembly_lookup[sample_name].keys():
			sys.exit("ERROR in assembly fofn" + assem_fof_name + " assembly file basenames must be unique")
		assembly_lookup[sample_name][abname] = assem_file
		assembly_types[sample_name][abname] = assem_type

for sample_name in reads.keys():
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".ccs.yak"))

alignments = {}
summaries = []
results_files = {}
for sample_name in assemblies.keys():
	haploid = 0
	alignments[sample_name] = {}
	results_files[sample_name] = []
	for assembly in assemblies[sample_name]:
		assem_bname = os.path.basename(assembly)
		alignments[sample_name][assem_bname] = []
		if assembly_types[sample_name][assem_bname] == "all":
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"yak.qv.out"))
			results_files[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"yak.qv.out"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.gaps.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.ev.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pchlst.ctg.bed"))
			results_files[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pchlst.ctg.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.primary.pchlst.ctg.bed"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"kat22.comp.stats"))
		if assembly_types[sample_name][assem_bname] == "all":
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"stats.txt"))
			results_files[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"stats.txt"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,assem_bname + ".full_table"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,assem_bname + ".duplication_summary"))
		results_files[sample_name].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,assem_bname + ".duplication_summary"))
		for reads_file in reads[sample_name]:
			read_bname = os.path.basename(reads_file)
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".bam"))
			if assembly_types[sample_name][assem_bname] == "all": 
				TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))
				alignments[sample_name][assem_bname].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))

	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".qc.summary.txt"))
	summaries.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".qc.summary.txt"))

TARGETS.append(os.path.join(config['output_base_dir'],"Results","Summary.txt"))

		
localrules: all, asset_step1, asset_step2, asset_step3, asset_step4, master_summary, summary

rule all:
	input: TARGETS

rule yak_hash:
	input:
		get_trimmed_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}" + ".ccs.yak")
	params:
		load = loadPreCmd(config["load"]["yak"])
	threads:
		32
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
		load = loadPreCmd(config["load"]["yak"])
	threads:
		32
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
		load = loadPreCmd(config["load"]["minimap2"])
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minimap2paf") * attempt
	shell:
		"{params.load} minimap2 -t {threads} -x asm20 {input[1]} {input[0]} > {output} || echo 'MINIMAP FAILED' > {output}"

rule minimap2bam:
	input:
		get_read_file,
		get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","{read_basename}" + ".bam")
	params:
		load = loadPreCmd(config["load"]["minimap2"],config["load"]["samtools"])
	threads:
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minimap2bam") * attempt
	shell:
		"{params.load} minimap2 -t {threads} -x asm20 {input[1]} {input[0]} > {output} || echo 'MINIMAP FAILED' > {output}"

rule asset_step1:
	input:
		assembly = get_assembly,
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.gaps.bed")
	params:
		load = loadPreCmd(config["load"]["asset"]),
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
		load = loadPreCmd(config["load"]["asset"]),
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
		load = loadPreCmd(config["load"]["asset"]),
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
		load = loadPreCmd(config["load"]["asset"]),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		1
	shell:
		"{params.load} cd {params.outdir}/{params.sample}/QC/{params.assem_bname} && pchlst -c {input[0]} {input[1]} > {output} || echo 'ASSET FAILED' > {output}"

rule asset_filt_primary:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.pchlst.ctg.bed")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","asset.primary.pchlst.ctg.bed")
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("asset_filt_primary") * attempt
	run:
		with open(output[0],"w") as outfile:
			with open(input[0], "r") as infile:
				for line in infile:	
					if "ASSET FAILED" in line:
						outfile.write("ASSET FAILED\n")
					elif "ptg" in line:
						outfile.write(line)
rule kat:
	input:
		reads = get_trimmed_reads,
		assembly = get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","kat22.comp.stats")
	params:
		load = loadPreCmd(config["load"]["kat"]),
		outdir = config['output_base_dir'],
		sample = get_sample,
		assem_bname = get_assem_bname
	threads:
		32
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("kat") * attempt
	shell:
		"{params.load} kat comp -t {threads} -m 22 -o {params.outdir}/{params.sample}/QC/{params.assem_bname}/kat22.comp '{input.reads}' {input.assembly} || echo 'KAT FAILED' > {output}"

rule asm_stats:
	input:
		get_assembly
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","stats.txt")
	params:
		load = loadPreCmd(config["load"]["bbmap"]),
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("asm_stats") * attempt
	shell:
		"{params.load} stats.sh {input} > {output}"

rule busco:
	input:
		get_assembly
	output:
		full_table = os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","{assembly_basename}.full_table"),
		short_sum = os.path.join(config['output_base_dir'],"{sample}","QC","{assembly_basename}","{assembly_basename}.duplication_summary")
	params:
		busco_db = lambda wc:get_busco_db(wc),
		outdir = config['output_base_dir'],
		load = loadPreCmd(config["load"]["busco"]),
		assem_bname = get_assem_bname
	threads:
		8
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco")
	shell:
		"{params.load} mkdir -p {params.outdir}/BUSCO_RUN/TMP && cd {params.outdir}/BUSCO_RUN && run_BUSCO.py -f -z -c {threads} -i {input} -l {params.busco_db} -o {params.assem_bname} -m genome -t {params.outdir}/BUSCO_RUN && cp {params.outdir}/BUSCO_RUN/run_{params.assem_bname}/short_summary_{params.assem_bname}.txt {output.short_sum} && cp {params.outdir}/BUSCO_RUN/run_{params.assem_bname}/full_table_{params.assem_bname}.tsv {output.full_table}"


rule summary:
	input:
		results = get_results_files
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.qc.summary.txt")
	params:
		sample_name = get_sample
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("summary") * attempt
	run:
		result_data = {}
		print("IN RUN")
		for filename in input.results:
			print("PICKED " + filename)
			result_basename = os.path.basename(filename)
			assembly_basename = os.path.basename(os.path.dirname(filename))
			if assembly_basename not in result_data.keys():
				result_data[assembly_basename] = {}
				result_data[assembly_basename]['qv'] = ""
				result_data[assembly_basename]['breaks'] = ""
				result_data[assembly_basename]['stats'] = ",,,,,"
			if result_basename == "yak.qv.out":
				print("FOUND YAK")
				result_data[assembly_basename]['qv'] = extract_qv(filename)
				print(result_data[assembly_basename]['qv'])
				print("DONE YAK")
			if result_basename == "asset.pchlst.ctg.bed":
				print("FOUND ASSET")
				result_data[assembly_basename]['breaks'] = str(sum(1 for _ in open(filename)))
				print(result_data[assembly_basename]['breaks'])
				print("DONE ASSETT")
			if result_basename == "stats.txt":
				print("FOUND STATS")
				result_data[assembly_basename]['stats'] = extract_stats(filename) 
				print(result_data[assembly_basename]['stats'])
				print("DONE STATS")
			if result_basename.endswith("duplication_summary"):
				print("FOUND BUSCO")
				result_data[assembly_basename]['busco'] = extract_busco(filename)
				print(result_data[assembly_basename]['busco'])
				print("DONE BUSCO")
				
		
		with open(output[0],"w") as outfile:
			outfile.write(",".join(["#Sample,Assembly,QV,NumContigs,AssemblySize,N50,N90,NumContigs>50kb,PercentAssembly>50kb,NumBreakpoints,BUSCO\n"]))
			outfile.write(",".join([params.sample_name,assembly_basename,result_data[assembly_basename]['qv'],result_data[assembly_basename]['stats'],result_data[assembly_basename]['breaks'],result_data[assembly_basename]['busco'] + "\n"]))

rule master_summary:
	input:
		summaries
	output:
		os.path.join(config['output_base_dir'],"Results","Summary.txt")
	params:
		outdir = config['output_base_dir'],
	threads:
		1
	shell:
		"mkdir -p {params.outdir}/Results && echo '#Sample,Assembly,QV,NumContigs,AssemblySize,N50,NumContigs>50kb,PercentAssembly>50kb,NumBreakpoints' > {output} && cat {input} | grep -v '#' >> {output}"
			
