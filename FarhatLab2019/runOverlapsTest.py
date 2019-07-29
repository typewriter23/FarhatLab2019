from countVCFOverlaps import countOverlaps
import os

# This script tests the injectivity of Omega^-1
    # I.E. tests the ability to call variants from RNA
    # When compared to DNA
    
    # Goes off of the **3SamplePaper**'s data
    
# dnaIDs = [SRR832991, SRR833024, SRR833121, SRR924700]
# RNAIDSets = 
inputFolder = "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs"

mappings = {"SRR832991": ["SRR6176429", "SRR6176430", "SRR6176437"],
            "SRR833024": ["SRR6176431", "SRR6176432", "SRR6176434"],
            "SRR833121": ["SRR6176433", "SRR6176435", "SRR6176436"],
            "SRR924700": ["SRR6176438", "SRR6176439", "SRR6176440"]}
            
def genFileName(runID):
    """Given a fileID, returns the name of its corresponding .vcf"""
    # TODO: Make this cleaner :P
    if runID in mappings.keys():
        return os.path.join(inputFolder, runID + "_1.vcf")
    else: 
        return os.path.join(inputFolder, runID + "_1.cleaned.vcf")
        
overlaps = {}
averages = {}
summary = {"both": 0,
            "dnaOnly": 0,
            "rnaOnly": 0,
            "total": 0}
for dna in mappings.keys():
    averages[dna] = {"both": 0,
            "dnaOnly": 0,
            "rnaOnly": 0,
            "total": 0}
    for rna in mappings[dna]:
        dnaFile = genFileName(dna)
        rnaFile = genFileName(rna)
        numOverlaps = countOverlaps(dnaFile, rnaFile) # Returns dictionary of overlaps
        overlaps[(dna, rna)] = numOverlaps
        # TO Done: Oh my god this is such bodged way to do things
            # I don't even know how to fix this though atm :/
            # Nvm I found a better way
        averages[dna]["both"] += numOverlaps["both"]
        averages[dna]["dnaOnly"] += numOverlaps["onlyX"]
        averages[dna]["rnaOnly"] += numOverlaps["onlyY"]
        averages[dna]["total"] += numOverlaps["all"]
    summary["both"] += averages[dna]["both"]
    summary["dnaOnly"] += averages[dna]["dnaOnly"]
    summary["rnaOnly"] += averages[dna]["rnaOnly"]
    summary["total"] += averages[dna]["total"]

print("rawer data = " + str(overlaps))
print("averages = " + str(averages))
print("summary = " + str(summary))
