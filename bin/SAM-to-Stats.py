# Generate TP, FP, TN, FN stats for hisat2 vs star test

import sys

SAMList = list(line.strip().split("\t") for line in open(sys.argv[1]))
simReadStats = list(line.strip().split("\t") for line in open(sys.argv[2]))
Output = open(sys.argv[3], "w")

goodreads = 0
badreads = 0
allreads = 0
for i in simReadStats:
    if i == simReadStats[0]:
        pass
    else:
        if "spikesequence" in i[0]:
            badreads = badreads + int(i[1])  #Summing all the reads that are 'spike-ins'. Should not map to DB
        else:
            goodreads = goodreads + int(i[1]) #Summing all reads that should map to DB
        allreads = allreads + int(i[1])

Total_TP = 0
Total_FP = 0
Total_TN = 0
Total_FN = 0
spiked = 0
readsMapped = 0

for samLine in SAMList:
    if samLine[0].startswith("@"): #if a header
        pass
    else:
        readsMapped = readsMapped + 1
        read = samLine[0].split("_simread")[0]
        mapped_to = samLine[2]
        if "spikesequence" in read: # spike-in read mapped to DB sequence (tRNAs and sno/miRNAs)
            Total_FP = Total_FP + 1
            spiked = spiked + 1
        elif mapped_to.startswith("chr"): # tRNA reads
            ### tRNAs have their own code chunk because I decided that any
            ### tRNA from a species mapping against any member of that 
            ### species is a true positive
            tRNA = mapped_to.split("-")[1]
            if tRNA in read:
                Total_TP = Total_TP + 1 
            else:
                Total_FP = Total_FP + 1
        else: # sno/miRNAs
            if mapped_to == read:
                Total_TP = Total_TP + 1
            else:
                Total_FP = Total_FP + 1
        
Total_TN = badreads - spiked # True negatives = the number of spike-in reads NOT mapped to the DB
Total_FN = goodreads - Total_TP # This number is all the positive (good) reads that were lost through mapping/collapse

print ("Max TP: ", goodreads, "Max TN: ", badreads)
print (Total_TP, Total_FP, Total_TN, Total_FN)

if int(Total_TP) != 0:
    Sensitivity = float(Total_TP)/(Total_TP + Total_FN)
if int(Total_TN) != 0:
    Specificity = float(Total_TN)/(Total_TN + Total_FP)

print ("Sensitivity: ", Sensitivity, "Specificity: ", Specificity)
Output.write("Total number of reads\t" + str(allreads) + "\nTrue positives:\t" + str(Total_TP) + "\nFalse positives:\t" + str(Total_FP) + "\nTrue negatives:\t" + str(Total_TN) + "\nFalse negatives:\t" + str(Total_FN) + "\nSensitivity:\t" + str(Sensitivity) + "\nSpecificity:\t" + str(Specificity) + "\n")	

Output.close()


