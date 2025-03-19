# suffixient-array

## Overview

`suffixient-array` is a software providing two main functionalities: 1) construction of a suffixient set of smallest cardinality $\chi$ for a text $T[1..n]$, 2) computation of the suffixient-array (sA) based index supporting locate one occurrence queries and MEMs finding as described in the paper "Suffixient Arrays: a New Efficient Suffix Array Compression Technique" ([go to the paper](https://arxiv.org/abs/2407.18753)).

## Download and Install

~~~~
git clone https://github.com/regindex/suffixient-array
cd suffixient-array
mkdir build
cd build
cmake ..
make
~~~~

## Usage
You can construct a smallest suffixient set of your input using the suffixient-set.py interface. This interface calls the C++ binaries that process the input text file and generate a suffixient set based on the chosen algorithm. Note that the input file must not contain the characters 0, 1, or 2, as they are reserved for internal use.
```
usage: suffixient-set.py [-h] [-a ALGORITHM] [-o OUTPUT] [-w WSIZE] [-p MOD] [-t THREADS] [-i] input

positional arguments:
  input                 input file name

options:
  -h, --help            show this help message and exit
  -a, --algorithm ALGORITHM
                        suffixient set construction algorithm (one-pass,linear,fm,PFP)
  -o, --output OUTPUT   output files basepath (def. input)
  -w, --wsize WSIZE     PFP: sliding window size (def. 10)
  -p, --mod MOD         PFP: hash modulus (def. 100)
  -t, --threads THREADS
                        PFP: number of helper threads (def. 0)
  -i, --no-invert       PFP: do not invert the text before running the algorithm
```

You can use the suffixient-array-index.py interface to build and query the suffixient-array index (sA-index). Like the previous tool, it requires an input file that must not contain the characters 0, 1, or 2. In addition, if the `opt-sA` index variant or the `rlz|bitpacked-text` random access text oracles  are selected the input must contain DNA characters (A,C,G,T) only.
```
usage: suffixient-array-index.py [-h] [-v INDEX_VARIANT] [-o ORACLE_VARIANT] [-b] [-e EXTRA_EF_SPACE] [-l] [-p PATTERN_FILE] input

positional arguments:
  input                 input files basepath

options:
  -h, --help            show this help message and exit
  -v, --index-variant INDEX_VARIANT
                        suffixient-array index variant: sA(def.)|opt-sA|PA
  -o, --oracle-variant ORACLE_VARIANT
                        random access text oracle variant: lz77(def.)|rlz|bitpacked-text
  -b, --build-index     build suffixient-array index
  -e, --extra-ef-space EXTRA_EF_SPACE
                        opt-sA: maximum allowed extra space (percentage) on top of the suffixient-array: Def. 30
  -l, --locate-one-occ  run locate one occurrence queries
  -p, --pattern-file PATTERN_FILE
                        file containing the patterns to locate (fasta format required)
```

### Run on Example Data

```console
// Construct a smallest suffixient set using the PFP compressed-space one-pass algorithm
python3 build/suffixient-set.py -a PFP data/yeast.fasta 

// Construct and query the baseline suffixient-array index (sA-index)
python3 build/suffixient-array-index.py --build-index data/yeast.fasta
python3 build/suffixient-array-index.py --locate-one-occ -p data/yeast_patterns.fasta data/yeast.fasta
```

### Datasets

You can download the datasets we used to evaluate the sA-index at the following link: https://github.com/regindex/suffixient-array/releases/download/datasets-bio-v1.0/biological.7z

### Requirements

The \texttt{suffixient-array} tool requires
* A Linux 64-bit operating system.
* A modern Python 3 release version 3.7 or higher.
* A modern C++11\14 compiler such as `g++` version 4.9 or higher.

## Reference and citation 

[1] Davide Cenzato, Francisco Olivares, Nicola Prezza: On Computing the Smallest Suffixient Set. SPIRE 2024: 73-87 ([go to the paper](https://doi.org/10.1007/978-3-031-72200-4_6))

[2] Davide Cenzato, Lore Depuydt, Travis Gagie, Sung-Hwan Kim, Giovanni Manzini, Francisco Olivares, Nicola Prezza: Suffixient Arrays: a New Efficient Suffix Array Compression Technique. CoRR abs/2407.18753 (2025) ([go to the paper](https://doi.org/10.48550/arXiv.2407.18753))

If you use \texttt{suffixient-array} in an academic setting, please cite this work as follows:
### suffixient-array
    @article{suffixient-array-2025,
      author       = {Davide Cenzato and
                      Lore Depuydt and
                      Travis Gagie and
                      Sung-Hwan Kim and
                      Giovanni Manzini and
                      Francisco Olivares and 
                      Nicola Prezza},
      title        = {Suffixient Arrays: a New Efficient Suffix Array Compression Technique},
      journal      = {CoRR},
      volume       = {abs/2407.18753},
      year         = {2025},
      doi          = {10.48550/ARXIV.2407.18753}
    }

## Funding

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe research and innovation programme, project REGINDEX, grant agreement No. 101039208.