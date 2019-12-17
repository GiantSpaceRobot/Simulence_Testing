from Bio import SeqIO
#from Bio.Seq import Seq
import sys

# Write all sequence names to file 
mySet = set()

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
out_file = open(sys.argv[2], "w")
for fasta in fasta_sequences:
    name = fasta.id
    if name.startswith("chr"):
        seqName = name.split("-")[1]
    else:
        seqName = name
    mySet.add(seqName) #Ensure no duplicates

for i in mySet:
    out_file.write(i + "\n")
out_file.close()

