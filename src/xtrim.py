#!/usr/bin/env python3

# Import libraries
import gzip, re, argparse

# Functions for entry processing
def control_entry(entry): 
    """Function ensures that the three first lines of each read are in the correct format. Returns True if they are valid, else returns False."""
    #Check that length of sequence and quality is the same
    if len(entry[1]) != len(entry[3]): 
        return False

    # Check the entry is in the correct format by checking line 1-3 
    if re.match(r"^@", entry[0]) is None: 
        return False

    if not re.match(r"^[AGCTXN]+$", entry[1]):
        return False

    if not re.match(r'^\+.*$', entry[2]):
        return False

    return True


def create_phred(N):
    """Creates a dictionary of the first 42 phred scores based on the number given as an input"""
    phred_dict = {}
    for i in range(0, 43):  
        phred_dict[chr(i + N)] = i

    return phred_dict


def convert_phred(qualentry, E = None):
    """Function takes the quality string as input, and encoding E, if given. Returns the string converted to Q scores in list, if string and E is valid.q"""
    Q_list = []

    # For phred quality encoding 33
    if E == 33:
        phred33 = create_phred(33)

        # Check the quality line of the entry and convert it to a list with converted quality scores
        for i in qualentry:
            if i in phred33:
                Q_list.append(int(phred33[i]))
            else: 
                return False
        return Q_list
    
    # For phred quality encoding 64
    elif E == 64: 
        phred64 = create_phred(64)

        # Check the quality line of the entry and convert it to a list with converted quality scores
        for i in qualentry:
            if i in phred64:
                Q_list.append(int(phred64[i]))
            else: 
                return False
        return Q_list
    
    # For the case that phred encoding is not given by user (it will be autodetected)
    elif E is None:  
        phred33 = create_phred(33)
        phred64 = create_phred(64)

        # Escape special characters and join the phred dicts into a regular expressioon
        REphred33 = '^' + '[' + ''.join(re.escape(char) for char in phred33.keys()) + ']' + '+$' 
        REphred64 = '[' + ''.join(re.escape(char) for char in phred64.keys()) + ']'
        
        # Check the quality line of the entry, and convert it to a list with converted quality scores
        if re.match(REphred33, qualentry): 
            for i in qualentry: 
                Q_list.append(int(phred33[i]))
            return Q_list
        
        # Check the quality line of the entry and convert it to a list with converted quality scores
        elif re.match(REphred64, entry):
            for i in qualentry: 
                Q_list.append(int(phred64[i]))
            return Q_list

        else: 
            return False
    
    # Controlling invalid phred encoding cases
    else: 
        return False


def trim(entry, Qlist, trimtype, thres3 = None, thres5 = None, mw = None):
    """Trims entries based on either number of bases N, or by quality threshold Q using moving window if window size is given."""
    
    # For trimming based on number of bases
    if trimtype == "N":
        if thres3 is not None:
            
            # Check that trim length does not exceed the length of the entry
            if thres3 >= len(entry[1]):
                return False
            
            # trim sequence and quality line in entry from 3' end
            entry[1] = entry[1][:-thres3]
            entry[3] = entry[3][:-thres3]

        if thres5 is not None: 
            # Check that trim length does not exceed the length of the entry
            if thres5 >= len(entry[1]):
                return False

            # trim sequence and quality line in entry from 5' end
            entry[1] = entry[1][thres5:]
            entry[3] = entry[3][thres5:]

        if thres3 is not None and thres5 is not None:
            # Check that trim length does not exceed the length of the entry
            if thres5+thres3 >= len(entry[1]):
                return False
        return entry
    
    # For trimming based on quality of bases
    elif trimtype == "Q":

        # if moving window size is not given, default is 1
        if mw is None:
           mw = 1
        
        # Check that moving window size does not exceed the length of the entry
        if mw > len(entry[1]):
            return False

        if thres5 is not None:
            i = 0

            # trim the entry form 5' by removing the first base in the moving window if mean quality of the moving window is lower than threshold
            while sum(Qlist[i:i+mw])/mw < thres5 and mw <= len(entry[1]):
                entry[1] = entry[1][1:]
                entry[3] = entry[3][1:]
                i += 1

        if thres3 is not None:
            j = 0

            # trim the entry form 3' by removing the last base in the moving window if mean quality of the moving window is lower than threshold
            while sum(Qlist[::-1][j:j+mw])/mw < thres3 and mw <= len(entry[1]):
                entry[1] = entry[1][:-1]
                entry[3] = entry[3][:-1]
                j += 1
        
        # if the entry is trimmed to be shorter than the length of the moving window, it will be discarded
        if len(entry[1]) < mw:
            return False

        return entry

    # check that valid trimtype is given
    else:
        raise ValueError("Trimtype should be N or Q.")

        
def postprocess(entry, trimQlist, length = None, qual = None, N = None):
    """Postprocessing of trimmed reads, checking for minimum length, minumum quality and maximum number of unknown bases."""
    flag = True

    # check read is longer than given threshold 
    if length is not None:
        if len(entry[1]) < length:
            flag = False
    
    # check mean trimmed read quality is higher than given threshold
    if qual is not None:
        if sum(trimQlist) / len(trimQlist) < qual: 
            flag = False

    # check number of N bases are below given threshold
    if N is not None:
        if entry[1].count('N') > N:
            flag = False
    
    # flag will be returned as true if the read passes all the requirements
    return flag


def write_log(kept, invalphred, invalentry, lowqual, overtrim, zipped, logfile): 
    """A log file is written, containing information such as the number of reads kept and discarded"""
    try:
        with open(logfile, 'w') as lf: 
            lf.write(f"The input file is a {'gzip' if zipped else 'fastq'} file \n")
            lf.write(f"Total number of entries: {kept + invalphred + invalentry + lowqual + overtrim} \n")
            lf.write(f"Number of trimmed entries: {kept} \n")
            lf.write(f"Number of entries with invalid phred quality: {invalphred} \n")
            lf.write(f"Number of invalid entries: {invalentry} \n")
            lf.write(f"Number of entries removed because of low quality after trimming (low mean quality, short length, N content): {lowqual} \n")
            lf.write(f"Number of entries removed because of invalid trimming parameters: {overtrim} \n")

    except FileNotFoundError:
        print(f"Input file not found.")
    except PermissionError:
        print(f"Permission denied for input file.")
    except Exception as e:
        print(f"An error occurred: {e}")


def main_process(entry, kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped, o):
    """Takes one entry at the time, and completes the process of checking and trimming, keeping track of how many entries are being kept, 
    and how many are being discarded, and write valid entries to outputfile"""

    # check entry is valid
    if control_entry(entry) is True:
        Qlist = convert_phred(entry[3], args.phred)
    
        # check that the enrty has a valid phred encoding
        if isinstance(Qlist, list):
            trimmedentry = trim(entry, Qlist, trimtype = args.trimtype, thres3 = args.thres3, thres5 = args.thres5, mw = args.movwin)

            # check that all trimming parameters are valid
            if isinstance(trimmedentry, list): 
                trimmedQlist = convert_phred(trimmedentry[3], args.phred)
                
                # check that the trimmed read satisfies the input requirements
                if postprocess(entry, trimmedQlist, args.minlen, args.minqual, args.maxN):
                    kept_reads += 1 #kept trimmed entry
                    
                    # write entry to output file
                    for line in trimmedentry:
                        o.write(line + "\n")
                else:
                    low_qual += 1 #trimmed entry did not satisfy input requirments
            else:
                overtrim += 1 #attempting to trim over entry's length
        else:
            inval_phred += 1 #phred is invalid
    else:
        inval_entry += 1 #entry is invalid

    return kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped


def readfile(filename, outfilename, args):
    """Opens files that are zipped or unzipped, and writes the file out zipped or unzipped accordingly."""

    inval_entry = inval_phred = overtrim = low_qual = kept_reads = 0

    # for gzip files
    try:
        with gzip.open(filename, 'rt') as f, gzip.open(outfilename, 'wt') as o:
            zipped = True

            while True:
                entry = [f.readline().strip() for _ in range(4)]  # Read four lines at a time
                if entry[0] == '':  # If the first line is empty, it means we've reached the end of the file and stop
                    break
                kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped = main_process(entry, kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped, o)
                write_log(kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped, args.log)

    except OSError:
        zipped = False

    # for fastq files
    if zipped == False:
        try: 
            with open(filename, 'r') as f, open(outfilename, 'w') as o:
                while True:
                    entry = [f.readline().strip() for _ in range(4)]  # Read four lines at a time
                    if entry[0] == '':  # If the first line is empty, it means we've reached the end of the file and stop
                        break
                    kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped = main_process(entry, kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped, o)
                    write_log(kept_reads, inval_phred, inval_entry, low_qual, overtrim, zipped, args.log)

        except FileNotFoundError:
            print(f"Input file not found.")
        except PermissionError:
            print(f"Permission denied for input file.")
        except Exception as e:
            print(f"An error occurred: {e}")



def main(infile, outfile, logfile, phred, trimtype, thres_3, thres_5, movwin, minlen, minqual, max_N):
    """main program function"""
    print("Input file:                                     ", infile)
    print("Output file:                                    ", outfile)
    print("Log file:                                       ", logfile)
    print("Phred-score:                                    ", phred)
    print("Trimming based on:                              ", trimtype)
    print("Threshold 3' end:                               ", thres_3)
    print("Threshold 5' end:                               ", thres_5)
    print("Moving window size:                             ", movwin)
    print("Minimum length after trimming:                  ", minlen)
    print("Minimum quality after trimming:                 ", minqual)
    print("Maximum number of unknown bases after trimming: ", max_N)

# Ascii Art 

xtrim_art = """  ░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░▒▓██████████████▓▒░  
  ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
   ░▒▓██████▓▒░   ░▒▓█▓▒░   ░▒▓███████▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░ """
                                                                     


print(xtrim_art)
print("              #########################################")
print("              #               XTrim v1.0              #")
print("              #########################################")
print(" ")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="XTrim: Readtrimmer for fastq files!")

    # Add command-line arguments
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file (Mandatory)")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file (Mandatory)")
    parser.add_argument("-lg", "--log", type=str, required=True, help="Log file (Mandatory)")
    parser.add_argument("-p", "--phred", type=int, required=False, help="Phred encoding (Optional)")
    parser.add_argument("-tt", "--trimtype", type=str, required=True, help="Trimming based on no. of bases(N), or quality(Q) (Mandatory)")
    parser.add_argument("-t3", "--thres3", type=int, required=False, help="No. of bases to trim / minimum quality threshold (Optional)")
    parser.add_argument("-t5", "--thres5", type=int, required=False, help="No. of bases to trim / minimum quality threshold (Optional)")
    parser.add_argument("-w", "--movwin", type=int, required=False, help="Size of moving window for trimming (Optional)")
    parser.add_argument("-l", "--minlen", type=int, required=False, help="Minimum length of read after trimming")
    parser.add_argument("-q", "--minqual", type=float, required=False, help="Minimum mean quality of read after trimming (Optional)")
    parser.add_argument("-N", "--maxN", type=int, required=False, help="Maximum number of unknow bases in read after trimming (Optional)")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with the provided inputs
    main(args.input, args.output, args.log, args.phred, args.trimtype, args.thres3, args.thres5, args.movwin, args.minlen, args.minqual, args.maxN)
    readfile(args.input, args.output, args)