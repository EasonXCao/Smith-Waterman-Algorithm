#!/usr/bin/python
__author__ = "Eason Cao"
__email__ = "first.last@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np
import pandas as pd

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

### runSW(inputFile, scoreFile, openGap, extGap)
### takes input sequence file, substitution matrix file, value for opening gap and extending gap as inputs;
### outputs the score matrix and the alignment result
def runSW(inputFile, scoreFile, openGap, extGap):
    
    # write the sequences as strings
    # write the substitution matrix in a dataframe score_df
    with open(inputFile) as inpf:
        seq = inpf.read().splitlines()
        seq1 = seq[0]
        seq2 = seq[1]
    score_df = pd.read_table(scoreFile, sep=' +', index_col=0, engine='python')
    
    ## initialize four matrices with same shape
    lc = len(seq1)
    lr = len(seq2)
    
    # S matrix represents the final score matrix
    S = np.zeros((lr+1, lc+1), dtype = int)

    # nogap matrix records all the scores without taking any gaps
    nogap = np.zeros((lr+1, lc+1), dtype = int)

    # gap1 matrix records all the scores by only making and extending gaps in sequence 1
    gap1 = np.zeros((lr+1, lc+1), dtype = int)

    # gap2 matrix records all the scores by only making and extending gaps in sequence 2
    gap2 = np.zeros((lr+1, lc+1), dtype = int)
    
    # path is a dictionary in which the key represents each position in S matrix and the value represents
    # the score of that position
    path = {}

    # insert a space in each sequence
    seq1m = ' '+seq1[:]
    seq2m = ' '+seq2[:]
    
    # iterate through all the positions in S matrix
    for row in range(1, lr+1):
        for col in range(1, lc+1):

            # calculate the score of the position for nogap, gap1 and gap2 based on the updated S matrix
            nogap[row, col] = S[row-1, col-1] + score_df.loc[seq2m[row], seq1m[col]]
            gap1[row, col] = max(np.add(S[0:row, col], (np.arange(row, 0, -1)-1)*extGap+openGap))
            gap2[row, col] = max(np.add(S[row, 0:col], (np.arange(col, 0, -1)-1)*extGap+openGap))

            # set the value in S matrix as the maximum of nogap, gap1, gap2 and 0
            S[row, col] = max([0, nogap[row, col], gap1[row, col], gap2[row,  col]])
            
    # record the path
    for i in range(0, lr+1):
        for j in range(0, lc+1):
            
            # the values in first row and first column of the S matrix come from nowhere
            if i==0 or j==0:
                path['[' + str(i) + ', ' + str(j) + ']'] = []
            else:
                path['[' + str(i) + ', ' + str(j) + ']'] = []
                
                # the remaining positions come from one of the postions of three matrices 
                if nogap[i, j] == S[i, j]:
                    path['[' + str(i) + ', ' + str(j) + ']'].append('[' + str(i-1) + ', ' + str(j-1) + ']')
                if gap1[i, j] == S[i, j]:
                    path['[' + str(i) + ', ' + str(j) + ']'].append('[' + str(i-1) + ', ' + str(j) + ']')
                if gap2[i, j] == S[i, j]:
                    path['[' + str(i) + ', ' + str(j) + ']'].append('[' + str(i) + ', ' + str(j-1) + ']')

    # write S matrix into a dataframe S_df
    S_df = pd.DataFrame(S, index = list(seq2m), columns = list(seq1m))
    
    # find the position of the largest number in the score matrix, which is the starting point of traceback
    end = np.argwhere(S == S.max())
    
    # if several positions have the same largest value, record them all as keys in path
    # value represents the value of the key
    # result represents the current key
    for i in end:
        key = str(list(i))
        value = path[key]
        result = [key]
        
    # output the result as a new txt file
    with open('output.txt', 'a') as output:
        output.write('-' * 11 + '\n')
        output.write('|Sequences|' + '\n')
        output.write('-' * 11 + '\n')
        output.write('sequence1' + '\n')
        output.write(seq1 + '\n')
        output.write('sequence2' + '\n')
        output.write(seq2 + '\n')

        output.write('-' * 14 + '\n')
        output.write('|Score Matrix|' + '\n')
        output.write('-' * 14 + '\n')
        output.write(S_df.to_csv(sep = '\t'))

        output.write('-' * 22 + '\n')
        output.write('|Best Local Alignment|' + '\n')
        output.write('-' * 22 + '\n')
        output.write('Alignment Score:' + str(int(S.max())) + '\n')
        output.write('Alignment Results:' + '\n')

        # call traceback(path, S, value, result, seq1m, seq2m)
        # which outputs the alignment result
        output.write(traceback(path, S, value, result, seq1m, seq2m))



### traceback(path, S, value, result, seq1, seq2)
### a recursive function that takes path, S, value, result, sequence1, sequence2
### outputs the alignment result
def traceback(path, S, value, result, seq1, seq2):
    
    # set the value as new key
    key = value[0]
    
    i = int((key.split(',')[0]).strip('['))
    j = int((key.split(',')[1]).strip(']'))
    
    if value != []:

        # when applying traceback, if there are several paths leading to the same position,
        # choosing the one with largest score
        if len(value) != 1:
            for l in range(0, len(value)):
                key_sub = value[l]
                i_sub = int((key_sub.split(',')[0]).strip('['))
                j_sub = int((key_sub.split(',')[1]).strip(']'))
                if S[i_sub, j_sub] > S[i, j]:
                    key = value[l]
        result.append(key)
        value = path[key]

    # the base case of the recursive function is the postion whose value is 0 
    if S[i, j] == 0:
        x = 0
        y = 0
        s1 = ''
        s2 = ''
        md = ''
        
        for n in range(len(result)-2, -1, -1):
            point = result[n]
            i = int((point.split(',')[0]).strip('['))
            j = int((point.split(',')[1]).strip(']'))
            
            # output the alignment result with given format
            if i == x:
                s1 += seq1[j]
                s2 += '-'
                md += ' '
            elif j == y:
                s1 += '-'
                s2 += seq2[i]
                md += ' '
            else:
                s1 += seq1[j]
                s2 += seq2[i]
                if seq1[j] == seq2[i]:
                    md += '|'
                else:
                    md += ' '
            if n == 0:
                s1 = s1 + ')' + seq1[j+1:]
                s2 = s2 + ')' + seq2[i+1:]
            elif n == len(result) - 2:
                pm1 = len(seq1[0:j])
                pm2 = len(seq2[0:i])
                pre_max = max(pm1, pm2)
                if pm1 == pre_max:
                    s1 = seq1[1:j] + '(' + s1
                    s2 = ' ' * (pm1 - pm2 - 1) + seq2[0:i] + '(' + s2
                    md = ' ' * pre_max + md
                elif pm2 == pre_max:
                    s1 = ' ' * (pm2 - pm1 - 1) + seq1[0:j] + '(' + s1
                    s2 = seq2[1:i] + '(' + s2
                    md = ' ' * pre_max + md
            x = i
            y = j
            
        return s1 + '\n' + md + '\n' + s2
    
    else:
        return traceback(path, S, value, result, seq1, seq2)

runSW(args.input, args.score, args.opengap, args.extgap)


## Usage: python hw1.py -i <input file> -s <score file>
## Example: python hw1.py -i input.txt -s blosum62.txt



