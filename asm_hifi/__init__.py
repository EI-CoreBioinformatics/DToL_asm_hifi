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
		assertFilenameValid(os.path.join(os.path.dirname(__file__), "etc", "asm_hifi.config.yaml"), "MISSING CONFIG ERROR: Failed to find default config")
	def create_run_config(self):
		with open(self.default_config_path,'r') as yamlfile:
			with open(self.run_config_path,'w') as run_config:
				for line in yamlfile:
					run_config.write(line)
				run_config.write("output_base_dir: " + self.args.outdir + "\n")
				run_config.write("sample_sheet: " + self.args.sample_sheet + "\n")
	def run(self):
		if os.path.exists(self.done_marker_path) and self.args.rerun:
			os.remove(self.done_marker_path)
		if not os.path.exists(self.done_marker_path):
			snakemake_success = run_snakemake(self.snakefile, self.snakemake_outdir, self.run_config_path, self.exe_env, dryrun=self.args.dryrun)
			if snakemake_success:
				Path(self.done_marker_path).touch()
				print("SUCCESS!")
			else:
				print("Snakemake terminated without success\n")
		else:
			print("Run is already complete")
	

class QCModule(AssemHiFi):
	def __init__(self,args):
		AssemHiFi.__init__(self,args)

class AssemblyModule(AssemHiFi):
	def __init__(self,args):
		AssemHiFi.__init__(self,args)
	def create_run_config(self):
		AssemHiFi.create_run_config(self)
		with open(self.run_config_path,'a') as run_config:
			if self.args.keepbin:
				run_config.write("HIFIASM_KEEP_BIN_FILES: True\n")
			else:
				run_config.write("HIFIASM_KEEP_BIN_FILES: False\n")

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
				run_config.write("assembly_sheet: " + self.args.sample_sheet + "\n")
