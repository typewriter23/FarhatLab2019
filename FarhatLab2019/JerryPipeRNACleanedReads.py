import os
import pandas as pd
import numpy as np
# from slurmpy import Slurm
# import vcf
import shutil
import ntpath
from JerryPipeToolbox import *

# JerryPipe, but with an added QC step that trims reads
    # Designed to process RNA sequencing data

###############################################################
########## Constants:##########################################
###############################################################
# envName = "JerryEnv" # Environment with all of the packages pre-loaded
    # # Packages included: {bwa, java, picard, pilon, samtools, ...}
# refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
    # # Absolute path to the H37RV genome
# scratchDir = "/n/scratch2/jy250/22SampleTest/"
# samFolder = "/n/scratch2/jy250/22SampleTest/samFiles/" 
# bamFolder = "/n/scratch2/jy250/22SampleTest/bamFiles/" 
# cleanedFastqFolder = "/n/scratch2/jy250/22SampleTest/cleanedFastqs/"
# inputFolder = "/n/scratch2/jy250/22SampleTest/"  # "/n/scratch2/jy250/test100/" # 
# vcfFolder = "/n/scratch2/jy250/22SampleTest/vcfFiles/"
# countsFolder = "/n/scratch2/jy250/22SampleTest/counts/"
# lineageFolder = "/n/scratch2/jy250/lineageCalls/"
# slurmFolder = "/n/scratch2/jy250/22SampleTest/slurmLogs/"

envName = "JerryEnv" # Environment with all of the packages pre-loaded
    # Packages included: {bwa, java, picard, pilon, samtools, ...}
refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
    # Absolute path to the H37RV genome
scratchDir = "/n/scratch2/jy250/1000SampleTest"
samFolder = "/n/scratch2/jy250/1000SampleTest/samFiles/" 
bamFolder = "/n/scratch2/jy250/1000SampleTest/bamFiles/" 
cleanedFastqFolder = "/n/scratch2/jy250/1000SampleTest/cleanedFastqs/"
inputFolder = "/n/scratch2/jy250/1000SampleTest"  # "/n/scratch2/jy250/test100/" # 
vcfFolder = "/n/scratch2/jy250/1000SampleTest"
slurmFolder = "/n/scratch2/jy250/1000SampleTest/slurmLogs/"
countsFolder = "/n/scratch2/jy250/1000SampleTest/counts/"
lineageFolder = "/n/scratch2/jy250/lineageCalls/"
# These functions will take names of files, and generate the coresponding 
        # name of the new file
        # e.g. abc.fasta ==> abc.sam
        
def genSamName(fastq): 
    """Function to systematically generate absolute paths to 
        sam file names from FASTQ file names"""
    return os.path.join(samFolder, os.path.splitext(fastq)[0] + ".sam")
    # return os.path.join(samFolder, ntpath.split(fastq)[1].replace(".fastq", ".sam"))
def genBamName(sam):
    return os.path.join(bamFolder, os.path.splitext(sam)[0] + ".bam")
def genVCFName(bam):
    return os.path.join(vcfFolder, os.path.splitext(bam)[0])
def genMetricsName(drbam):
    return os.path.join(vcfFolder, os.path.splitext(drbam)[0] + ".metrics")
def genBaseName(fileName):
    """Removes the '_1' and '.fastq' extentsions from a filename"""
    return fileName.split("_")[0].split(".")[0]
def genSampleID(path):
    """Generates the sampleID from the filename given
        Code partially taken from StackOverflow:
        https://stackoverflow.com/questions/8384737/extract-file-name-from-path-no-matter-what-the-os-path-format"""
    head, tail = ntpath.split(path)
    result =  tail or ntpath.basename(head)
    return result.split(".")[0] # Gets just the sample name, cleans out the ".cleaned.[EXT]"
def genSLURMName(path):
    return os.path.join(slurmFolder, os.path.splitext(path)[0] + ".slurm")
# def genSamName(fastq): 
    # """Function to systematically generate absolute paths to 
        # sam file names from FASTQ file names"""
    # return os.path.join(samFolder, os.path.splitext(fastq)[0] + ".sam")
    # # return os.path.join(samFolder, ntpath.split(fastq)[1].replace(".fastq", ".sam"))
# def genBamName(sam):
    # return os.path.join(bamFolder, genSampleID(sam) + ".bam")
# def genVCFName(bam):
    # return os.path.join(vcfFolder, genSampleID(bam))
# def genMetricsName(drbam):
    # return os.path.join(vcfFolder, genSampleID(drbam) + ".metrics")
# def genBaseName(fileName):
    # """Removes the '_1' and '.fastq' extentsions from a filename"""
    # return fileName.split("_")[0].split(".")[0]
# def genSampleID(path):
    # """Generates the sampleID from the filename given
        # Code partially taken from StackOverflow:
        # https://stackoverflow.com/questions/8384737/extract-file-name-from-path-no-matter-what-the-os-path-format"""
    # head, tail = ntpath.split(path)
    # result =  tail or ntpath.basename(head)
    # return result.split(".")[0] # Gets just the sample name, cleans out the ".cleaned.[EXT]"
###############################################################
###############################################################

def genPairedFileName(fileName):
    """Creates the name of the paired end file, assuming that one is given a file that ends with "_1.fastq":
        XYZ_1.fastq => XYZ_2.fastq"""
    # TODO: Find a cleaner way to do this
    return fileName.replace("_1.fastq", "_2.fastq")
    
listOfFastqs = [] # ["SRR6176440_1.fastq", "SRR6176434_1.fastq"]
for filename in os.listdir(inputFolder):
    if filename.endswith("_1.fastq"): 
        pairedReads = [filename, genPairedFileName(filename)]
        listOfFastqs.append(pairedReads)
        
tag = "tag"
output_dir = "/n/scratch2/jy250/output/"
scratch_dir = "/n/scratch2/jy250/"
   
# def setUpSLURM(commands_list):
    # Nothing
    # Adds the heading of a slurm script, does nothing for now
def runQC(commands_list, fastq, fastqPaired = None, minQual = 10):
    """ Quality control steps, including: {prinseq read trimming}"""
    # TODO: Add QC steps
    # TODO: Add Kraken
    return trimReads(commands_list, fastq, fastqPaired, minQual)
    
def trimReads(commands_list, fastq, fastqPaired = None, minQual = 5):
    """Function to trim reads using prinseq
        
        Takes in a (pair of) FASTQ file(s) and appends to commands_list 
         the commands needed to trim each read of the files to a minimum phred 
         score of MINQUAL
         
        Outputs the name of the cleaned FASTQ file or files"""
         
    leftMinQual = minQual
    rightMinQual = minQual
    outputFile = os.path.join(inputFolder, genBaseName(fastq) + ".cleaned") # Give it an absolute path # TODO: HARDCODED! BAD
        #TODO: Hardcoded to write into the input folder => change later!!!
    if fastqPaired is not None:
        trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                -out_good {outputFileName} -fastq {fastq1} -fastq2 {fastqPaired}".format(leftMinQual = leftMinQual, 
                rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq, fastqPaired = fastqPaired)
        commands_list.append(trimCommand)
    else:
        trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                -out_good {outputFileName} -fastq {fastq1}".format(leftMinQual = leftMinQual, 
                rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq)
        commands_list.append(trimCommand)
    return genCleanedOutputName(outputFile, paired = True)
    
def Launch_JerryPipe(fastqFile, pairedFile = None, minQ = 10):
    '''
    This script is based on Max Marin's JankyPipe. Repurposed to work with Conda and RNAseq data.
    
    Converts a single (pair of) FASTQ file into {a VCF file, a COUNTS.TSV file, and a LINEAGE.TXT file} 
    
    "This script launches a job to call variants for the input fastq files against H37Rv
    using a number of packages. The important output (VCF, lineage info files, quality report)
    is stored in the output directory while the intermediary files (SAMs, trimmed fastqs, BAM, etc)
    are stored in a scratch directory."
    
    - Fundamental components:
        - Generate a file (a SLURM script) with the following commands:
            -- take a FASTQ file, and convert it to SAM via (bwa)
            -- take the SAM file and convert it to BAM via (samtools)
            -- *** take the SAM file and convert it to VCF via (Pilon) ***
            
        - For all of the files, submit the job via SBATCH command
    - Store all commands in a list
    - Submit commands to slurm
    '''
    commands_list = []
    
    # setUpSLURM(commands_list) # Probably delete later, set up the requirements for the slurm script to run
    prepConda(commands_list)
    fastqFileCleaned = [fastqFile, pairedFile] # runQC(commands_list, fastqFile, pairedFile, minQual = minQ)
    #Skipping cleaning steps b/c we're only working with cleaned data
    
    # We just want to keep track of the names
    #TODO: Add functionality for paired end reads
    samFile = convertToSAM(commands_list, refGenome, fastqFileCleaned[0], fastqFileCleaned[1]) # Adds commands to commandList to convert 
                                                      # the named fastq file to the named sam file
    bamFile = convertToBAM(commands_list, samFile)                 # Adds commands to commandList to convert 
                                                      # the named sam file to the named bam file
    cleanBamFile = cleanBams(commands_list, bamFile)
    vcfFile = convertToVCF(commands_list, cleanBamFile, refGen = refGenome, outputDir = "/") 
    countReads(commands_list, cleanBamFile)
    
    # "\n".join(commands_list)
    # print(commands_list)
    slurmOutput = os.path.join(slurmFolder, genSLURMName(fastqFile))
    submitSlurmScript(commands_list, slurmOutput)
   
for fastqPair in listOfFastqs:
    # Handling for non-paired end data
    try:
        fileName = os.path.join(inputFolder, fastqPair[0])
        pairedFile = os.path.join(inputFolder, fastqPair[1])
        Launch_JerryPipe(fileName, pairedFile)
    except IOError: 
        fileName = os.path.join(inputFolder, fastqPair[0])
        pairedFile = os.path.join(inputFolder, fastqPair[1])
        Launch_JerryPipe(fileName)
# Launch_JerryPipe("/n/scratch2/jy250/fastqFiles/SRR6176430_1.fastq")
