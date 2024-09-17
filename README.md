XTrim: Readtrimmer for fastq files
Introduction
XTrim is a Python-based tool designed for trimming fastq files. It provides functionality to trim and evaluate reads based on various criteria such as quality scores, length and ambiguous base content.
Installation
XTrim requires Python3 to run. To install, simply download the provided python script and ensure that you have the necessary libraries installed. You can install the required libraries using pip.
Usage
The tool is to be run from command line as such:
python xtrim.py −−input <inputfile> −−output <outputfile> −−logfile <logfile> −−trimtype <Q/N> [−−phred <33/64>] [−−thres3 <threshold 3’>]
[−−thres5 <threshold 5’>] [−−movwin <moving window size>][−−minlen <minimum length>] [−−minqual <minimum mean quality>][−−maxN <maximum N content>]
or
python xtrim.py −i <inputfile> −o <outputfile> −lg <logfile> −tt <Q/N> [−p <33/64>] [−t3 <threshold 3’>] [−t5 <threshold 5’>] [−w <moving window size>]
[−l <minimum length>] [−q <minimum mean quality>] [−N <maximum N content>]
The program needs to be run from the src directory.
