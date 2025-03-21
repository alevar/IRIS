# IRIS: Detection and Validation Of Chimeric Reads.

[![PyPI version](https://badge.fury.io/py/iris-av.svg)](https://pypi.org/project/iris-av/)
[![GitHub Downloads](https://img.shields.io/github/downloads/alevar/iris/total.svg)](https://github.com/alevar/IRIS/releases/latest)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Introduction

IRIS is a method designed to detect and validate chimeric junction from multi-genome alignments. The method constructs a DP alignment matrix
from two separate alignments to infer precise breakpoint. The two-pass algorithm is implemented to refine consistency of breakpoint inference.
The method is designed to take advantage of anntoations of either or, ideally, both genomes involved in the chimeric event by penalizing and prioritizing events at known junctions.

## Publications

Coming soon...

## Documentation

## Installation

### Via PyPI

The easiest way to install IRIS is through PyPI:

```bash
$ pip install iris-av
$ iris --help
```

To uninstall SNAPPER:

```bash
$ pip uninstall iris-av
```

### Building from source

To build from source, clone the git repository:

```bash
$ git clone https://github.com/alevar/iris.git --recursive
$ cd iris
$ pip install -r requirements.txt
$ pip install .
```

### Requirements

| Requirement | Details |
| ----------- | ------- |
| Language support | Python ≥ 3.6 |
| Dependencies | - |

## Getting started

IRIS expects BLAST alignments to be provided in the following format:

```bash
blastn \
  -db blast_database \
  -query query_fasta \
  -out output.blastn6 \
  -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

### Usage

```bash
iris [-h] -i1 INPUT1 -i2 INPUT2 -a1 ANNOTATION1 -a2 ANNOTATION2 -o OUTPUT [--two_pass] [-g] [--chim-genome]
                   [-g1 GENOME1] [-g2 GENOME2] [-max_dist MAX_DIST] [-max_weight MAX_WEIGHT]
                   [-full_weight FULL_WEIGHT] [-half_weight HALF_WEIGHT] [--overhang OVERHANG]
```

### Options

| Option | Description |
| ------ | ----------- |
| `-i1, --input1` | Path to the file containing BLAST mapping of reads to genome 1. Alignment is expected to have the following format: -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" |
| `-i2, --input2` | Path to the file containing BLAST mapping of reads to genome 2. Alignment is expected to have the following format: -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" |
| `-a1, --annotation1` | Path to the file containing GTF/GFF annotation for genome 1. |
| `-a2, --annotation2` | Path to the file containing GTF/GFF annotation for genome 2. |
| `-g1, --genome1` | Path to the file containing genome 1 FASTA sequence. |
| `-g2, --genome2` | Path to the file containing genome 2 FASTA sequence. |
| `--two-pass` | Flag enables the 2-pass mode. Breakpoints from the first pass will be used to bias DP trace towards consensus sites. |
| `--group` | If enabled, will output a file with breakpoints groupped by position. |
| `-o, --output` | Path to the output file. |
| `--chim-genome` | (Requires -group). If enabled, will generate a fasta file with chimeric genome sequences, stitching together the two genomes at the breakpoints. |
| `--max-dist` | Maximum distance between breakpoints of the two segments. Default: 5. |
| `---max-weight` | Maximum weight of a breakpoint when biasing the 2nd pass. Default: 5. |
| `--full-weight` | Weight of a breakpoint that matches donor and acceptor. Default: 5. |
| `--half-weight` | Weight of a breakpoint that matches either donor or acceptor. Default: 3. |
| `--overhang` | Number of bases to include in the chimeric genome overhang. Default: 1000. |

### Help Options

| Option | Description |
| ------ | ----------- |
| `-h, --help` | Prints help message. |

## Example Data

Sample datasets are provided in the "example" directory to test and get familiar with SNAPPER.

The included example can be run with the following command from the root directory of the repository:

```bash
iris --input1 ./examples/AY69_E4p5_LTA/host.blastn.6 --input2 ./examples/AY69_E4p5_LTA/path.blastn.6 --annotation1 ./examples/AY69_E4p5_LTA/host.gtf --annotation2 ./examples/csess.1.0.0.known.gtf --output ./examples/AY69_E4p5_LTA/ris --genome1 ./examples/AY69_E4p5_LTA/host.fa --genome2 ./examples/SIV239.fa --chim-genome --two-pass --group
```