#!/bin/bash -e
rm -rf dist/
rm -rf asm_hifi.egg-info
rm -rf build
python setup.py bdist_wheel
pip install -I --prefix=/ei/workarea/group-ga/Projects/CB-GENANNO-488_DToL_Assembly_polishing_pipeline/Analysis/TestBuild/x86_64 -U dist/*whl
