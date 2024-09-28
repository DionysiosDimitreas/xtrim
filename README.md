# XTrim: Readtrimmer for FASTQ Files

**XTrim** is a Python-based tool designed for trimming FASTQ files. It trims sequencing reads based on quality scores, read length, and ambiguous base content.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
  - [Mandatory Arguments](#mandatory-arguments)
  - [Optional Arguments](#optional-arguments)
- [Example Usage](#example-usage)
- [Contributing](#contributing)

## Introduction
XTrim is a flexible tool that helps in preprocessing FASTQ files by trimming and evaluating sequencing reads based on:
- Quality scores (Phred)
- Read length
- Ambiguous base content (`N` bases)

The tool ensures cleaner, more reliable datasets for downstream analysis.

## Installation
XTrim requires **Python 3**. Install the required dependencies using the following command:

```bash
pip install -r requirements.txt

−−maxN/−N: Maximum number of unknown bases in sequence after trimming\
\
Example Usage:\
python xtrim.py −i input.fastq −o output.fastq −lg logfile.log −tt N −t3 6 −t5 8−l 50−q 32−N 7
```
## Usage

XTrim is run via the command line and accepts various arguments for trimming sequencing reads. You can use either a long or short form for the arguments.

Long format
```bash
python xtrim.py --input <inputfile> --output <outputfile> --logfile <logfile> --trimtype <Q/N> \
[--phred <33/64>] [--thres3 <3’end threshold>] [--thres5 <5’end threshold>] \
[--movwin <moving window size>] [--minlen <min length>] [--minqual <min quality>] [--maxN <max N content>]
```

Short format
```bash
python xtrim.py -i <inputfile> -o <outputfile> -lg <logfile> -tt <Q/N> \
[-p <33/64>] [-t3 <3’end threshold>] [-t5 <5’end threshold>] \
[-w <moving window size>] [-l <min length>] [-q <min quality>] [-N <max N content>]
```
Mandatory arguments:\
−−input/−i: input file, gz or txtfile (same format as output)\
−−output/−o: output file , gz or txtfile (same format as input)\
−−logfile/−lg: log file\
−−trimtype/−tt: Q or N, trimming type, defining whether to trim based on
length or quality Optional arguments:\
−−phred/−p : 33 or 64, phred encoding\
−−thres3/−t3: Threshold 3’−end\
−−thres5/−t5: Threshold 5’−end\
−−movwin/−w: Size of moving window for quality trimming\
−−minlen/−l: Minimum length of read after trimming\
−−minqual/−q: Minimum mean quality of read after trimming\
