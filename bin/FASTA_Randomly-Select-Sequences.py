from Bio import SeqIO
#from Bio.Seq import Seq
import sys
import random

# Randomly select sequences from a file 

# When the 10% chance for adding tRNA sequences to spike-in file is true, 
# this code checks whether a member of this tRNA species has been added to 
# spike-in file first. If it was already added, it does not add it again.
# This greatly simplifies the calculation of TP, FP, TN, FN at the stats step later.


fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
spike_file = open(sys.argv[3] + "_spike.fa", "w")
notspike_file = open(sys.argv[3] + "_not-spike.fa", "w")
percentage_of_sequences = int(sys.argv[2]) # The % of sequences you want to randomly select from the input
tRNAset = set()
tRNAlist = list()

for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    if "chr" in name:
        tRNAspecies = name.split("-")[1]
        tRNAlist.append(tRNAspecies)

spikeSet = set()

### Choose 10 tRNA species at random. These will be in the spiked dataset
for i in range(10):  
    tRNAspecies = random.choice(tRNAlist)
    spikeSet.add(tRNAspecies)

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
for fasta in fasta_sequences:
    randomNum = random.randint(1,100)
    name, sequence = fasta.id, fasta.seq
    #tRNAspecies = "AStringThatDoesNotOccurNormally" 
    if "chr" in name:
        tRNAspecies = name.split("-")[1]
        if tRNAspecies in spikeSet:
            spike_file.write(">spikesequence_" + name + "\n" + str(sequence) + "\n")
        else:
            notspike_file.write(">" + name + "\n" + str(sequence) + "\n")
    else:
        if randomNum < percentage_of_sequences:
            spike_file.write(">spikesequence_" + name + "\n" + str(sequence) + "\n")
        else:
            notspike_file.write(">" + name + "\n" + str(sequence) + "\n")

spike_file.close()
notspike_file.close()
