# Collapses SAM files

import sys

SAMList = list(line.strip().split("\t") for line in open(sys.argv[1]))
Output = open(sys.argv[2], "w")

myDict = dict()

for samLine in SAMList:
    line = "\t".join(samLine)
    if samLine[0].startswith("@"): #if a header
        Output.write(line + "\n")
    #elif samLine[0].startswith("chr"):
    else:
        read = str(samLine[0])#.split("_simread")[0])
        if read in myDict.keys():
            myDict[read].append(line)
        else:
            myList = [line]
            myDict[read] = (myList) 
    #else:
    #    Output.write(line + "\n")    

for k,v in myDict.iteritems():
    if k.startswith("chr"):
        feature = k.split("-")[1].split("_")[0]
    else:
        feature = k.split("_")[0]
    database_matches = list()
    for i in v:
        match = i.split("\t")[2]
        database_matches.append(match)
    #print (k, 
    if all(feature in x for x in database_matches) == True:  # if every value in the dict contains the feature species
        Output.write(v[0] + "\n")     

Output.close()
