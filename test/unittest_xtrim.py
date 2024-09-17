#!/usr/bin/env python3

import pytest
import sys

sys.path.append('src')

from xtrim import *

def test_createphred1():
    assert create_phred(33) == {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41, 'K': 42}, "Check that function creates the correct phred-dict based on given phred encoding"

def test_controlentry1():
    assert control_entry(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./']) == True, "Check that a valid entry is correctly translated to a list of quality scores"

def test_controlentry2():
    assert control_entry(['Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./']) == False, "Check that the function returns false in case the header is in wrong format"

def test_controlentry3():
    assert control_entry(['@Header1', 'ACCTGAACGNAAXT', '--', '!"#$%&()*+,-./']) == False, "Check that the function returns false in case the third line does not consist of a '+' sign"

def test_controlentry4():
    assert control_entry(['@Header1', 'ACCDDUYQCGTATW', '+', '!"#$%&()*+,-./']) == False, "Check that the function returns false in case the sequence contains invalid characters"

def test_controlentry5():
    assert control_entry(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()-.']) == False, "Check that the function returns false in case the sequence has a different length than the quality line in the entry"

def test_convertphred1():
    assert convert_phred('!"#$%&()*+,-./', 33) == [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Check that function correctly converts phred33 encoding correctly"

def test_convertphred2():
    assert convert_phred('!"#$%&()*+,-./') == [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Check that function correctly detects and converts phred33 encoding correctly"

def test_convertphred3():
    assert convert_phred('!"#$%&()*+,-./', 64) == False, "Check that function returns False when incorrect phred encoding has been provided"

def test_convertphred4():
    assert convert_phred('!"#$%&()*+,-./', 99) == False, "Check that function returns False for invalid phred encoding numbers"

def test_convertphred5():
    assert convert_phred('@ABCDEFGHIJKLMN', 64) == [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], "Check that function converts phred64 encoding correctly"

def test_trim1():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "N", thres3 = 2, thres5 = 3) == ['@Header1', 'TGAACGNAA', '+', '$%&()*+,-'] , "Check that function trims correctly when trimming with a specified number of reads"

def test_trim2():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "N", thres3 = 2) == ['@Header1', 'ACCTGAACGNAA', '+', '!"#$%&()*+,-'], "Check that function trims correctly when trimming with a specified number of reads"

def test_trim3():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "N", thres5 = 20) == False, "Check that the function returns False when the trim size is larger than the length of the read"
 
def test_trim4():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "N", thres3 = 10, thres5 = 10) == False, "Check that the function returns False when the trim size is larger than the length of the read"
    
def test_trim5():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Q", thres5 = 5) == ['@Header1', 'AACGNAAXT', '+', '&()*+,-./'], "Check that the function trims correctly when given quality treshold for only 5'"

def test_trim6():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Q", thres5 = 3, mw = 3) == ['@Header1', 'CTGAACGNAAXT', '+', '#$%&()*+,-./'], "Check that function trims correctly using moving window when trimming by quality"

def test_trim7():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Q", thres5 = 3, mw = 20) == False, "Check that function returns False for moving window size larger than the length"

def test_trim8():
    assert trim(['@Header1', 'ACCTGAACGNAAXTGG', '+', '!"#$%&()*+,-./!#'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14,0,2], "Q", thres3 = 2, thres5 = 3, mw = 2) == ['@Header1', 'TGAACGNAAXTG', '+', '$%&()*+,-./!'], "Check that function trims correctly from 5' and 3' end with moving window"

def test_trim9():
    assert trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], "Q", thres5 = 20, mw = 3) == False, "Check that function returns False for quality thresholds higher than the sequence quality"

def test_trim10():
    # Check that trimtype is either N or Q
    with pytest.raises(ValueError):
       trim(['@Header1', 'ACCTGAACGNAAXT', '+', '!"#$%&()*+,-./'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14], 90, thres5 = 20) 

def test_postproc1():
    assert postprocess(['@Header1', 'ACCTGAACGNAAXTGG', '+', '!"#$%&()*+,-./!#'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14,0,2], length = 10, qual = 3, N=3) == True, "Check that the function returns True when the read satisfies all specified requirements"
    
def test_postproc2():
    assert postprocess(['@Header1', 'ACCTGAACGNAAXTGG', '+', '!"#$%&()*+,-./!#'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14,0,2], length = 20, qual = 3, N=3) == False, "Check that the function discards reads that are shorter than specified"

def test_postproc3():
    assert postprocess(['@Header1', 'ACCTGAACGNAAXTGG', '+', '!"#$%&()*+,-./!#'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14,0,2], length = 10, qual = 20, N=3) == False, "Check that the function discards reads with lower mean quality than specified by user"

def test_postproc4():
    assert postprocess(['@Header1', 'ACCTNNNNGNAAXTGG', '+', '!"#$%&()*+,-./!#'], [0,1,2,3,4,5,7,8,9,10,11,12,13,14,0,2], length = 10, qual = 3, N=2) == False, "Check that the function discards reads with to many N in the sequence"

# Random text file
# 1 of the reads have diff seq length than quality length
# Entries with diff phred score values
# Unit testing:
# 1) Is it enough to only test the individual functions?
# 2) Do we need to test with many files and stuff like that and how do we hand that in?
# 3) When unit testing our .py file, we get an error that could not recognize 'args'. ???
# 4) dow we need to do unit testing on writing file functions? THAT WOULD BE VERY TIME CONSUMING