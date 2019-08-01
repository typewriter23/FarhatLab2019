# Script to find non-paired end info
from JerryPipeToolbox import *

def find(list1, list2):
    return list(set(list1) - set(list2))

# inputFolder = 
# list1 = []
# for filename in os.listdir(inputFolder):
    # if filename.endswith("_1.fastq"): 
        # list1.append(genSampleID(filename)
        
list1 = []
list2 = []

idListLocation1 = "../pairedEndRunIDs.txt"
with open(idListLocation1) as idFile:
    for line in idFile:
        list1.append(genSampleID(line))
        
idListLocation2 = "../pairedEndCompleted.txt"
with open(idListLocation2) as idFile:
    for line in idFile:
        list2.append(genSampleID(line))

remaining = find(list1, list2)
with open('remainingPairedEnd.txt', 'w') as f:
    for line in remaining:
        f.write(line)
