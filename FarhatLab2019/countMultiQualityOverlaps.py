from JerrPipeRNAPairedEnds import *

# Script that runs JerryPipe with varying miniumum phred scores for 
    # read trimming
envName = "JerryEnv" # Environment with all of the packages pre-loaded
refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
scratchDir = "/n/scratch2/jy250/"
samFolder = "/n/scratch2/jy250/dnaSamFiles/" 
bamFolder = "/n/scratch2/jy250/dnaBamFiles/" 
inputFolder = "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs" # "/n/scratch2/jy250/22SampleInputs/"
vcfFolder = "/n/scratch2/jy250/dnaVcfFiles/"
countsFolder = "/n/scratch2/jy250/counts"

minQuality = 1 # phred score for lowest quality run
maxQuality = 30 # phred score for highest quality run
listOfQualities = list(range(minQuality, maxQuality))
for qual in listOfQualities:
    runJerryPipe(qual)

def runJerryPipe(qual): 
    for fastqPair in listOfFastqs:
        fileName = os.path.join(inputFolder, fastqPair[0])
        pairedFile = os.path.join(inputFolder, fastqPair[1])
        Launch_JerryPipe(fileName, minQ = qual)
