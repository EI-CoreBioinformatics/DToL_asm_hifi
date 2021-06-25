# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

name="asm_hifi"
version = "0.1"

if sys.version_info.major != 3:
	raise EnvironmentError("""asm_hifi is a python module that requires python3, and is not compatible with python2. Also, it is now 2019 and support for 2.x will cease soon.""")

setup(
	name=name,
	version=version,
	description=description,
	long_description=long_description,
	author="Gareth Linsmith",
	author_email="gareth.linsmith@earlham.ac.uk",
	license="MIT",
	classifiers=["Development Status :: 4 - Beta","Topic :: Scientific Engineering :: Bio/Informatics","License :: OSI Approved :: MIT License","Operating System :: POSIX :: Linux",'Programming Language :: Python :: 3.4',"Programming Language :: Python :: 3.5","Programming Language :: Python :: 3.6"],
	zip_safe=False,
	keywords="genome assembly QA bioinformatics",
	packages=find_packages(exclude=["test"]),
	install_requires=["snakemake>=4.4.0","drmaa","matplotlib"],
	entry_points={"console_scripts": ["asm_hifi=asm_hifi.__main__:main"]},
	include_package_data=True,
	data_files=[('asm_hifi/etc',['asm_hifi/etc/asm_hifi.config.yaml','asm_hifi/etc/hpc_config.json']),('asm_hifi/zzz',['asm_hifi/zzz/QCModule.smk.py','asm_hifi/zzz/AssemblyModule.smk.py'])]
)

