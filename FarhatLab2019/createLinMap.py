from JerryPipeToolbox import callLineage, genSampleID
import csv
import os
# Script to run through all .LINEAGE.TXT files
    # and create a dictionary that maps RunID
    # to lineage called
    
inputFolder = "/n/scratch2/jy250/lineageCalls"

mappings = {}
for filename in os.listdir(inputFolder):
    if filename.endswith(".lineage.txt"):
        file = os.path.join(inputFolder, filename)
        ls = []
        with open(file, newline='') as csvfile:
            csvReader = csv.reader(csvfile, delimiter=',')
            ls = [lin for lin in csvReader]
        sampleID = genSampleID(filename)
        mappings[sampleID] = ls
        if (len(ls)!= 0):
            print(str(sampleID) + "\t" + str(mappings[sampleID][0][0]))
    # print(mappings)