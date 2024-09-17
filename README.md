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

Mandatory arguments:
−−input/−i: input file, gz or txtfile (same format as output) −−output/−o: output f i l e , gz or txtfile (same format as input) −−logfile/−lg: log f i l e
−−trimtype/−tt: Q or N, trimming type, defining whether to trim based on
length or quality Optional arguments:
−−phred/−p : 33 or 64, phred encoding
−−thres3/−t3: Threshold 3’−end
−−thres5/−t5: Threshold 5’−end
−−movwin/−w: Size of moving window f o r quality trimming
−−minlen/−l: Minimum length of read after trimming
−−minqual/−q: Minimum mean quality of read after trimming
−−maxN/−N: Maximum number of unknown bases in sequence after trimming
Example Usage:
python xtrim.py −i input.fastq −o output.fastq −lg logfile.log −tt N −t3 6 −t5 8−l 50−q 32−N 7
