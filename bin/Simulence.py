#!/usr/bin/env python

#Copyright 2019 Paul Donovan
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this 
#software and associated documentation files (the "Software"), to deal in the Software 
#without restriction, including without limitation the rights to use, copy, modify, merge, 
#publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
#persons to whom the Software is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or 
#substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
#PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
#OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
#OTHER DEALINGS IN THE SOFTWARE.

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"

#Import libraries
from Bio import SeqIO
from random import *
import sys
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Format (e.g. FASTA/FASTQ)')
parser.add_argument('FASTA input')
parser.add_argument('Read length (e.g "50")')
parser.add_argument('Fold coverage (e.g "100")')
parser.add_argument('FASTA/FASTQ output')
if len(sys.argv[1:]) == 0:
    """
    Simulence generates FASTA/FASTQ with realistic coverage distribution
    """
    parser.print_help()
    parser.exit()
args = parser.parse_args()

#Define parameters
readLength = int(sys.argv[3])
foldCov = float(sys.argv[4])
binaryList = [0,1]
ten_list = [0,1,0,0,0,0,0,0,0,0]
nucleotides = "AGTC"

#Define functions
def repeat_to_length(string_to_expand, length):
   return (string_to_expand * ((length/len(string_to_expand))+1))[:length]

def chunkIt(seq, num):
    # From: https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
    avg = len(seq) / float(num)
    reads = []
    last = 0.0
    while last < len(seq):
        reads.append(seq[int(last):int(last + avg)])
        last += avg
    return reads

def mutator(my_read):
    #Introduce mutations/INDELs into reads (error profiles from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4787001/)
    readString = ""
    for nucleotide in my_read:

        ### Mutation:
        mutateFloat = random()
        if mutateFloat < 0.0035:
            new_nuc = nucleotide
            while str(new_nuc) == str(nucleotide):
                new_nuc = choice(nucleotides)
            nucleotide = new_nuc
        else:
            pass
        
        ### INDELs:
        indelFloat = random()
        if indelFloat < 3.95e-6:   # equivalent to 3.95 x 10^-6
            if (choice(binaryList) == 1): # Flip a coin for heads/tails. If heads, insertion. If tails, deletion.   
                insertion = choice(nucleotides)
                nucleotide = nucleotide + insertion
            else: 
                nucleotide = ""
        
        ### Write (possibly) modified read to string
        readString = readString + nucleotide
    return readString
 
def read_generator(seq, readLen, foldCoverage):
    read_list = list()
    seq1 = seq
    seqLen = len(seq1)
    readsFor1X = float(seqLen)/float(readLen) #Number of reads required to get 1X coverage
    for n in range(int(foldCoverage)):   # for 1... in 100 read coverage
        flag = 1
        my_seq = seq1
        seqLen = len(seq1)
        newReadLen = int(readLen)
        random_choice1 = choice(binaryList) # Choose which side of the sequence to simulate this batch of reads from
        random_choice2 = choice(ten_list) # Randomly choose 1 out of 10
        for i in range(int(readsFor1X)):     ### Randomly select sequence from either end of given seq
            if random_choice1 == 1: # 50:50. Simulate reads from left-hand side of sequence
                new_read = my_seq[0:newReadLen] # Start new read from left of input sequence (randomly add 1 to both half of the time for additional randomness)
                my_seq = my_seq[newReadLen:] # Shorten seq1 to reflect removal of newly simulated read (sampling without replacement)
            else:                         # 50:50. Sim reads from right-hand side of sequence
                new_read = my_seq[-newReadLen:] # Start new read from right of input sequence
                my_seq = my_seq[:-newReadLen] # Shorten seq1 to reflect removal of newly simulated read (sampling without replacement)
            new_read = mutator(new_read)
            read_list.append(new_read)
    return read_list

fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
out_file = open(sys.argv[5], "w")
log_file = open(str(sys.argv[5]) + ".log", "w") 
### Write header for log file
log_file.write("Sequence\tNumber of reads generated\n")

for fasta in fasta_sequences:
    readLength = sys.argv[3] #reset this variable with every iteration
    name, sequence, description = fasta.id, str(fasta.seq), str(fasta.description)
    reads = read_generator(sequence, readLength, foldCov)
    count = 1
    newList = (str(name), str(len(reads)), "\n")
    log = "\t".join(newList)
    log_file.write(log)
    for new_seq in reads:
        new_name = name + "_simread" + str(count)
        readLength = len(new_seq)
        if sys.argv[1] == "FASTA":
            out_file.write(">" + new_name + "\n" + new_seq + "\n") #FASTA output
        elif sys.argv[1] == "FASTQ":
            quality_score = repeat_to_length("I", readLength)  # Using Q-score of 40 (I) for simulated FASTQ file
            out_file.write("@" + new_name + "\n" + new_seq + "\n+\n" + quality_score + "\n") #FASTQ output
        count = count + 1
out_file.close()
log_file.close()
