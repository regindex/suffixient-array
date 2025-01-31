# Compute the suffixient set index

### Overview

Software implementing the suffixient set index.

### Install

~~~~
git clone https://github.com/regindex/suffixient-index
cd suffixient-index
mkdir build
cd build
cmake ..
make
~~~~

### Run

Construct the suffixient index for --input and compute all MEMs for the reads in --mems: 

~~~~
python3 mems_suffixient.py --input ../data/paper_example.txt --mems ../data/paper_reads.fasta
~~~~

### Funding

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe research and innovation programme, project REGINDEX, grant agreement No. 101039208.