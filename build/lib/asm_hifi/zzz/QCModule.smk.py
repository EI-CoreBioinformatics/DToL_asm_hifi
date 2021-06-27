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

	return num_contigs,assem_size,L50,L90,num_over50kb,perc_over50kb

HPC_CONFIG = HpcConfig(DEFAULT_HPC_CONFIG_FILE)
TARGETS = []
assemblies = {}
reads = {}
reads_lookup = {}
assembly_lookup = {}
adapter_logs = []
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
		adapter_logs.append(reads_file + ".log")
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
summaries = []
#adapter_logs = []
for sample_name in assemblies.keys():
	haploid = 0
	alignments[sample_name] = {}
	for assembly in assemblies[sample_name]:
		assem_bname = os.path.basename(assembly)
		alignments[sample_name][assem_bname] = []
		if "L0" in assem_bname:
			haploid = 1
		if "all.contigs" in assem_bname:
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"yak.qv.out"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.gaps.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.ev.pb.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.pchlst.ctg.bed"))
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"asset.primary.pchlst.ctg.bed"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"kat22.comp.stats"))
		if "p_ctg" in assem_bname:
			TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,"stats.txt"))
		for reads_file in reads[sample_name]:
			#TARGETS.append(reads_file + ".log.parse")
			#adapter_logs.append(reads_file + ".log")
			read_bname = os.path.basename(reads_file)
			if "all.contigs" in assem_bname: 
				TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))
				alignments[sample_name][assem_bname].append(os.path.join(config['output_base_dir'],sample_name,"QC",assem_bname,read_bname + ".paf"))

	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".collate.adapter.logs"))
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".adapter.summary.txt"))

	if haploid == 1:
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".haploid.qc.summary.txt"))
		summaries.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".haploid.qc.summary.txt"))
	else:
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".qc.summary.txt"))
		summaries.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".qc.summary.txt"))
		TARGETS.append(os.path.join(config['output_base_dir'],sample_name,"QC",sample_name + ".qc.report.html"))

TARGETS.append(os.path.join(config['output_base_dir'],"Results","Summary.txt"))

print(adapter_logs)
		
localrules: all, asset_step1, asset_step2, asset_step3, asset_step4, master_summary, parse_adapter_summary, sum_adapter, summary, summary2html

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
		load = loadPreCmd(config["load"]["yak"])
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
		load = loadPreCmd(config["load"]["minimap2"])
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
		16
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("kat") * attempt
	shell:
		"{params.load} kat comp -t {threads} -m 22 -o {params.outdir}/{params.sample}/QC/{params.assem_bname}/kat22.comp '{input.reads}' {input.assembly} || echo 'KAT FAILED' > {output}"

rule parse_adapter_summary:
	input:
		adapter_logs
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.collate.adapter.logs")
	threads:
		1
	shell:
		"cat {input} | grep \"Total reads processed\|Reads with adapters\" | sed 's/ //g' | sed 's/,//g' | cut -f 1 -d '(' | sed 's/:/\t/' > {output}"

rule sum_adapter:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.collate.adapter.logs")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.adapter.summary.txt")
	threads:
		1
	run:
		total_reads = 0
		total_adapters = 0
		with open(output[0],"w") as outfile:
			with open(input[0], "r") as infile:
				for line in infile:
					fields = line.rstrip().split("\t")
					if "Totalreadsprocessed" in line:
						num_reads = int(fields[1])
						total_reads = total_reads + num_reads
					if "Readswithadapters" in line:
						num_adapters = int(fields[1])
						total_adapters = total_adapters + num_adapters
			outfile.write(str(total_reads) + "\t" + str(total_adapters) + "\n")
			
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

rule summary_haploid:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L0.all.contigs.fasta","yak.qv.out"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L0.all.contigs.fasta","asset.primary.pchlst.ctg.bed"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L0.p_ctg.fasta","stats.txt"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.adapter.summary.txt")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.haploid.qc.summary.txt")
	params:
		sample_name = get_sample
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("summary_haploid") * attempt
	run:
		with open(input[0], "r") as yak:
			for line in yak:
				if line.startswith("QV"):
					fields = line.rstrip().split("\t")
					qvL0 = fields[0]
		breaksL0 = sum(1 for _ in open(input[1]))
		num_contigs_L0, assem_size_L0, L50_L0, L90_L0, num_over50kb_L0, perc_over50kb_L0 = extract_stats(input[2])

		with open(input[3], "r") as adapt:
			for line in adapt:
				fields = line.rstrip().split("\t")
				num_reads = fields[0]
				reads_adapt = fields[1]

		with open(output[0],"w") as outfile:
			outfile.write(",".join(["#Sample,NumReads,NumReadsAdapter,PurgeLevel,QV,NumContigs,AssemblySizeize,N50,N90,NumContigs>50kb,PercentAssembly>50kb,NumBreakpoints\n"]))
			outfile.write(",".join([params.sample_name,num_reads,reads_adapt,"L0",qvL0,num_contigs_L0, assem_size_L0, L50_L0, L90_L0, num_over50kb_L0, perc_over50kb_L0,breaksL0 + "\n"]))

rule summary:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L1.all.contigs.fasta","yak.qv.out"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L1.all.contigs.fasta","asset.primary.pchlst.ctg.bed"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L1.p_ctg.fasta","stats.txt"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L2.all.contigs.fasta","yak.qv.out"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L2.all.contigs.fasta","asset.primary.pchlst.ctg.bed"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.L2.p_ctg.fasta","stats.txt"),
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.adapter.summary.txt")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.qc.summary.txt")
	params:
		sample_name = get_sample
	threads:
		1
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("summary") * attempt
	run:
		with open(input[0], "r") as yak:
			for line in yak:
				if line.startswith("QV"):
					fields = line.rstrip().split("\t")
					qvL1 = fields[0]
				elif "YAK FAILED" in line:
					qvL1 = "YAK FAILED"
		with open(input[3], "r") as yak:
			for line in yak:
				if line.startswith("QV"):
					fields = line.rstrip().split("\t")
					qvL2 = fields[0]
				elif "YAK FAILED" in line:
					qvL2 = "YAK FAILED"

		breaksL1 = str(sum(1 for _ in open(input[1])))
		breaksL2 = str(sum(1 for _ in open(input[4])))

		num_contigs_L1, assem_size_L1, L50_L1, L90_L1, num_over50kb_L1, perc_over50kb_L1 = extract_stats(input[2])
		num_contigs_L2, assem_size_L2, L50_L2, L90_L2, num_over50kb_L2, perc_over50kb_L2 = extract_stats(input[5])
		
		with open(input[6], "r") as adapt:
			for line in adapt:
				fields = line.rstrip().split("\t")
				num_reads = fields[0]
				reads_adapt = fields[1]
		
		with open(output[0],"w") as outfile:
			outfile.write(",".join(["#Sample,NumReads,NumReadsAdapter,PurgeLevel,QV,NumContigs,AssemblySizeize,N50,N90,NumContigs>50kb,PercentAssembly>50kb,NumBreakpoints\n"]))
			outfile.write(",".join([params.sample_name,num_reads,reads_adapt,"L1",qvL1,num_contigs_L1, assem_size_L1, L50_L1, L90_L1, num_over50kb_L1, perc_over50kb_L1,breaksL1 + "\n"]))
			outfile.write(",".join([params.sample_name,num_reads,reads_adapt,"L2",qvL2,num_contigs_L2, assem_size_L2, L50_L2, L90_L2, num_over50kb_L2, perc_over50kb_L2,breaksL2 + "\n"]))

rule summary2html:
	input:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.qc.summary.txt")
	output:
		os.path.join(config['output_base_dir'],"{sample}","QC","{sample}.qc.report.html")
	threads:
		1
	run:
		a = pd.read_csv(input[0]) 
		a.to_html(output[0]) 
		html_file = a.to_html()

rule master_summary:
	input:
		summaries
	output:
		os.path.join(config['output_base_dir'],"Results","Summary.txt")
	threads:
		1
	shell:
		"echo '#Sample,PurgeLevel,QV,NumContigs,AssemblySizeize,N50,NumContigs>50kb,PercentAssembly>50kb,NumBreakpoints' > {output} && cat {input} | grep -v '#' >> {output}"
			
