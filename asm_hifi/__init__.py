import pkg_resources
import sys
import os
import yaml
import errno
import collections
import shutil
from pathlib import Path
from abc import ABC, abstractmethod
from .snakemake_helper import *
from .utils import *
import json
from os.path import *

__title__ = "asm_hifi"
__author__ = "Gareth Linsmith (garethlin)"
__email__ = "gareth.linsmith@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2019 Earlham Institute"
__version__ = pkg_resources.require("asm_hifi")[0].version

#REQUIRED BY SNAKEMAKE HELPER ExecutionEnviroment constructor
DEFAULT_HPC_CONFIG = assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "hpc_config.json"),"MISSING CONFIG ERROR: Failed to find HPC config")

class Organism:
        def __init__(self, name, kingdom, closest_reference=None, transcriptome=None):
                self.name = name
                self.kingdom = kingdom
                self.closest_reference = closest_reference
                self.transcriptome = transcriptom
class SeqLibrary:
	def __init__(self, name):
		self.name = name

class CssCell(SeqLibrary):
	def __init__(self,movname,filepath):
		self.movname = movname
		self.filepath = filepat

class AssemblyProject:
        def __init__(self,name,organism,libs = {},assemblies = {}):
                self.name = name
                self.organism = organism
                self.libs = libs
                self.primary_lib = None
                self.assemblies = assemblies

class SampleSheet:
	def __init__(self,filename):
		self.filepath = filepath
		self.samples = {}
		self.read()
	def read(self):
		infile = open(self.filename)
		for line in infile:
			if line.startswith('#'):
				continue
			fields = line.rstrip().split(",")
			if len(fields) == 3:
				organism_id = fields[0]
				movname = fields[1]
				readsfile = fields[2]
				kingdom = None
				closest_reference = None
				transcriptome = None
			#if len(fields) == 4:
			#add code for closest ref.. etc..
			else:
				sys.exit("FATAL ERROR: Invalid number of fields in sample sheet. " + str(len(fields)) + " fields detected in " + filename)
			

			if organism_id in self.samples.keys():
				mycell = CssCell(movname,readsfile)
				self.samples[organism_id].libs.append(mycell)
			else:	
				myorg = Organism(organism_id, kingdom, closest_reference, transcriptome)
				myproj = AssemblyProject(organism_id,myorg)
				mycell = CssCell(movname,readsfile)
				myproj.libs.append(mycell)
				self.samples[organism_id] = myproj

class AssemblySheet():
	def __init__(self,filename):
		self.filename = filename
		self.assemblies = {}
		self.read()
	def read(self):
		infile = open(self.filename)
		for line in infile:
			if line.startswith('#'):
				continue
			fields = line.rstrip().split(",")
			if len(fields) != 2:
				sys.exit("Invalid number of fields in assemblies file " + self.filename + "\n" + "Correct format is sample_name,assembly_path")
			else:
				sample_name = fields[0]
				assembly_path = fields[1]
				if " " in sample_name or "\t" in sample_name:
					sys.exit("Error in assemblies file " + self.filename + " spaces not alowed in sample names")
				if " " in assembly_path  or "\t" in assembly_path:
					sys.exit("Error in assemblies file " + self.filename + " spaces not alowed in assembly paths")
				if sample_name in self.assemblies.keys():
					sys.exit("Error: Duplicate sample name in assemblies file " + self.filename)
				else:
					if assembly_path in self.assemblies.values():
						sys.exit("Error: Duplicate assembly path in assemblies file " + self.filename)
					else:
						self.assemblies[sample_name] = assembly_path
	def get_assemblies():
		return self.assemblies

class AssemHiFi:
	def __init__(self,args):
		self.args = args
		self.snakefile = assertFilenameValid(os.path.join(os.path.dirname(__file__), "zzz", self.__class__.__name__ + ".smk.py"), "MISSING SNAKEFILE ERROR")
		self.default_config_path = assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "asm_hifi.config.yaml"), "MISSING CONFIG ERROR: Failed to find default config")
		self.hpc_logdir = createDirIfNotExist(os.path.join(self.args.outdir, "hpc_logs"))
		self.snakemake_outdir = createDirIfNotExist(os.path.join(self.args.outdir, "snakemake"))
		self.run_config_path = os.path.join(self.snakemake_outdir,self.__class__.__name__ + ".run.config")
		self.done_marker_path = os.path.join(self.snakemake_outdir,self.__class__.__name__ + ".done")
		self.read_default_config()
		self.create_run_config()
		self.exe_env = ExecutionEnvironment(self.args, NOW, job_suffix="assem", log_dir=os.path.join(self.hpc_logdir))
	def read_default_config(self):
		pass
	def create_run_config(self):
		pass
	def run(self):
		if not os.path.exists(self.done_marker_path):
			snakemake_success = run_snakemake(self.snakefile, self.snakemake_outdir, self.run_config_path, self.exe_env, dryrun=self.args.dryrun)
			if snakemake_success:
				#Path(self.done_marker_path).touch()
				print("SUCCESS!")
			else:
				print("Snakemake terminated without success\n")
		else:
			print("Run is already complete")
	

class PolishingModule(AssemHiFi):
	def __init__(self,args):
		AssemHiFi.__init__(self,args)
	def read_default_config(self):
		assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "asm_hifi.config.yaml"), "MISSING CONFIG ERROR: Failed to find default config")
	def create_run_config(self):
		with open(self.default_config_path,'r') as yamlfile:
			with open(self.run_config_path,'w') as run_config:
				for line in yamlfile:
					run_config.write(line)
				run_config.write("output_base_dir: " + self.args.outdir + "\n")
				run_config.write("assembly_sheet: " + self.args.assembly_sheet + "\n")

class QCModule(AssemHiFi):
	def __init__(self,args):
		AssemHiFi.__init__(self,args)
	def read_default_config(self):
		assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "asm_hifi.config.yaml"), "MISSING CONFIG ERROR: Failed to find default config")
	def create_run_config(self):
		with open(self.default_config_path,'r') as yamlfile:
			with open(self.run_config_path,'w') as run_config:
				for line in yamlfile:
					run_config.write(line)
				run_config.write("output_base_dir: " + self.args.outdir + "\n")
				run_config.write("sample_sheet: " + self.args.sample_sheet + "\n")

class AssemblyModule(AssemHiFi):
	def __init__(self,args):
		AssemHiFi.__init__(self,args)
	def read_default_config(self):
		assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "asm_hifi.config.yaml"), "MISSING CONFIG ERROR: Failed to find default config")
	def create_run_config(self):
		with open(self.default_config_path,'r') as yamlfile:
			with open(self.run_config_path,'w') as run_config:
				for line in yamlfile:
					run_config.write(line)
				run_config.write("output_base_dir: " + self.args.outdir + "\n")
				run_config.write("sample_sheet: " + self.args.sample_sheet + "\n")
				run_config.write("HIFIASM_KEEP_BIN_FILES: False\n")
