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

### Usage

```bash
iris [-h] -s SAM -r REFERENCE [-o OUTPUT] [--qry_intron_match_score QRY_INTRON_MATCH_SCORE] 
        [--trg_pos_match_score TRG_POS_MATCH_SCORE] [--trg_pos_mismatch_score TRG_POS_MISMATCH_SCORE]
```

### Options

| Option | Description |
| ------ | ----------- |
| `-s, --sam` | Path to the SAM/BAM alignment file. Read names in the alignment are expected to match corresponding transcript_id in the reference annotationPath to the query GTF/GFF annotation file. |
| `-r, --reference` | Path to the reference annotation. transcript_id field is expected to match read names in the sam alignment. |
| `-o, --output` | Path to the output SAM/BAM file. |
| `--qry_intron_match_score` | Score for matching query introns. |
| `--trg_pos_match_score` | Score for matching target positions. |
| `--trg_pos_mismatch_score` | Score for mismatching target positions. |

### Help Options

| Option | Description |
| ------ | ----------- |
| `-h, --help` | Prints help message. |

## Example Data

Sample datasets are provided in the "example" directory to test and get familiar with SNAPPER.

The included example can be run with the following command from the root directory of the repository:

```bash
iris --sam example/example.gtf --reference example/example.gtf --output example/output.sam
```