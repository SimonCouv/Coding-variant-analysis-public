#!/usr/bin/env bash
source activate dtp-rotation2
conda install git pip

# optional: clone latest versions of dfply and pybiomart
# note this may break the script if breaking changes are introduced in these packages
#git clone https://github.com/kieferk/dfply.git
#git clone https://github.com/jrderuiter/pybiomart.git

# install dfply
cd dfply
pip install -e .

# install pybiomart
cd ../pybiomart
pip install -e .



