import csv
import sys
import yaml
import argparse
import os
import pathlib
import subprocess

from collections import OrderedDict

from . import __version__
from .__init__ import *
from .snakemake_helper import *

def add_default_options(parser):
	common_group = parser.add_argument_group("asm_hifi options")
	common_group.add_argument("--outdir", "-o", type=str, default="asm_hifi")
	common_group.add_argument("--prefix", type=str, default="asm_hifi")
	common_group.add_argument("--dryrun", help="dry run",action="store_true")
	make_exeenv_arg_group(parser, default_hpc_config_file=DEFAULT_HPC_CONFIG, allow_mode_selection=False, silent=True)

def add_assembly_parser(subparsers):
	assembly_parser = subparsers.add_parser(
		"Assemble",
		help="",
		description=""
	)
	assembly_parser.add_argument("sample_sheet", type=str)
	add_default_options(assembly_parser)
	assembly_parser.set_defaults(runmode="assemble")

def add_qc_parser(subparsers):
	qc_parser = subparsers.add_parser(
		"QC",
		help="",
		description=""
	)
	qc_parser.add_argument("sample_sheet", type=str)
	qc_parser.add_argument("--assembly_path", type=str)
	add_default_options(qc_parser)
	qc_parser.set_defaults(runmode="qc")

def main():
	print("Starting EI asm_hifi V " + __version__)
	if len(sys.argv) == 1:
		sys.argv.append("-h")
	ap = argparse.ArgumentParser(prog="asm_hifi", description="The Earlham Institute HIFI Data Genome Assembly Pipeline (asm_hifi)")

	subparsers = ap.add_subparsers(
		help=""
	)

	add_assembly_parser(subparsers)
	add_qc_parser(subparsers)
	
	args = ap.parse_args()
	if args.runmode == "assemble":
		module = AssemblyModule(args)
	elif args.runmode == "qc":
		module = QCModule(args)
	else:
		sys.exit("Module not recognised")
	module.run()

if __name__ == "__main__":
	main()
