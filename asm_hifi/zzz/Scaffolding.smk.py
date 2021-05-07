import os.path

def get_sample(wc):
	return wc.sample

def get_assembly_path(wc):
	return assemblies[wc.sample]

def get_assembly_basename(wc):
	return os.path.basename(assemblies[wc.sample])

def get_reads(wc);
	return reads[wc.sample]

localrules: all, mkdir_and_copy_asm, mkref


TARGETS = []
assemblies = {}
reads = {}

infile = open(config['assembly_sheet'])
for line in infile:
	fields = line.rstrip().split(",")
	sample_name = fields[0]
	assembly_path = fields[1]
	readsdir = fields[2]
	assemblies[sample_name] = assembly_path
	reads[sample_name] = readsdir
	TARGETS.append(os.path.join(config['output_base_dir'],sample_name,sample_name + ".longranger.done"))
	#make this longranger.done to pull longranger rule

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
	shell:
		"mkdir -p {params.outdir}/{params.sample} && cp {input} {output}"

rule longranger:
	input:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa"),
		get_reads
	output:
		os.path.join(config['output_base_dir'],"{sample}","{sample}.longranger.done")
	params:
		load = "source longranger-2.1.2",
		outdir = config['output_base_dir'],
		sample = get_sample,
		reads = get_reads
	shell:
		"{params.load} && cd {params.outdir}/{params.sample} && longranger basic --id {params.sample} --fastqs {params.reads} && touch {output}"

#rule testpipe:
#	input:
#		os.path.join(config['output_base_dir'],"{sample}","{sample}.fa")
#	output:
#		os.path.join(config['output_base_dir'],"{sample}","{sample}.done")
#	shell:
#		"touch {output}"
