import os
import glob
from countVCFOverlaps import countOverlaps, countOverlapsVCFGZVCF
# Script that searches the database of 10,000 
    # MTB genomes for the closest match for a given
    # RNAseq VCF
    
rnaVCF = "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs/SRR6176432_1.cleaned.vcf"
inputFolder = "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data"

# def getAllVCFs(inputFolder):
    # """Function that returns an iterator over all vcf files
        # in the directory"""
    # for folder in os.listdir(inputFolder):
        # try:    
            # vcfFile = os.path.join(inputFolder, "pilon/")
            # for potentialFile in os.listdir(vcfFile):
                # if(potentialFile.endswith(".vcf")):
                    # yield potentialFile
        # except:
            # continue
            # for files in os.listdir(str(folder) + "/pilon/")

def getAllVCFs(inputFolder):
    """Function that returns an iterator over all vcf files
        in the directory:
        /n/data1/hms/dbmi/farhat/rollingDB/genomic_data 
        (Designed specifically with this O2 directory in mind)"""
    listOfFiles = []
    numErrors = 0
    for folder in os.listdir(inputFolder): #Loop through all folders
        # print(folder)
        try: # using a try-xcept block to avoid errors with other files 
            # not following this structure
            # TODO: make this cleaner
            vcfLoc = os.path.join(inputFolder, os.path.join(folder, "pilon/"))
            # print(vcfLoc)
            for potentialFile in os.listdir(vcfLoc):
                # print(potentialFile)
                if(potentialFile.endswith(".vcf.gz")):
                    listOfFiles.append(os.path.join(vcfLoc, potentialFile))
        except:
            # print("error at " + folder)
            numErrors += 1
    print(numErrors)
    return listOfFiles
def innerProd(vcfResults):
    """Function that maps results of countOverlaps (i.e. {both = , onlyX =, ...})
        into a real number
        
        Created for more flexibility (i.e. in case I want to change my distance metric"""
    both = vcfResults.get("both", 0)
    onlyX = vcfResults.get("onlyX", 0)
    onlyY = vcfResults.get("onlyY", 0)
    
    # return (both - onlyX - onlyY)/(both + onlyX + onlyY) 
    return (both)/(both + onlyX + onlyY)
        # Distance heuristic = # correct variants / total 
               # => % correct variants called
    
def compare(vcfX, gzvcfY, innerProd):
    """'Inner product' function for vcf files
        
        Given 2 vcf files ({vcfX, vcfY}) and an inner product function
        (innerProd), calculates the distance between them
        using distance = # of SNPs different"""
    regions = countOverlapsVCFGZVCF(vcf = vcfX, gzvcf = gzvcfY)
        # countOverlaps(vcfX, vcfY)
    return innerProd(regions)
    # TODO: figure out a distance metric
    # return regions["both"]/(regions["onlyX"] + regions["onlyY"] + regions["both"])

def compareToRef(ref, innerProdFun = innerProd):
    """Returns a 1-variable 
            function that returns compare([insertX], ref, innerProd)"""
    def compareFun(x):
        return compare(ref, x, innerProdFun)
    return compareFun
        
def getMin(listOfVCFs, compareFun, numMins = 1):
    # """Returns the numMin keys with smallest values in the list"""
    """Returns the smallest values in the list w.r.t. 
        a reference and a defined inner product fuction"""
    return min(listOfVCFs, key = compareFun)

print(getAllVCFs(inputFolder))
db = list(getAllVCFs(inputFolder)) # List of all vcfs in folder
minimumVCF = getMin(db, compareToRef(rnaVCF, innerProd), numMins = 1)
print(countOverlaps(rnaVCF, minimumVCF))

