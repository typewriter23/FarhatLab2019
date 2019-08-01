import os
import pandas as pd
import numpy as np
# from slurmpy import Slurm
# import vcf
import shutil
import ntpath

# Functions that are called by JerryPipe
    
# These 3 functions will take names of files, and generate the coresponding 
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
    return genBaseName(result.split(".")[0]) # Gets just the sample name, cleans out the ".cleaned.[EXT]"
def genLinFileName(vcfFile):
    """Generates the file location to which to write the lineage information
        "XXX.cleaned_1.vcf" => XXX.lineage.txt"""
    return os.path.join(lineageFolder, genSampleID(vcfFile) + ".lineage.txt")
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
def genReadcountsName(bamFile):
    return os.path.join(countsFolder,  genSampleID(bamFile) + ".counts.tsv")
    #   genBaseName(vcfFile) + ".counts.txt"

###############################################################
###############################################################
# envName = "JerryEnv" # Environment with all of the packages pre-loaded
# refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
# scratchDir = "/n/scratch2/jy250/"
# samFolder = "/n/scratch2/jy250/dnaSamFiles/" 
# bamFolder = "/n/scratch2/jy250/dnaBamFiles/" 
# inputFolder = "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs" # "/n/scratch2/jy250/22SampleInputs/"
# vcfFolder = "/n/scratch2/jy250/dnaVcfFiles/"
# countsFolder = "/n/scratch2/jy250/counts"

# envName = "JerryEnv" # Environment with all of the packages pre-loaded
    # # Packages included: {bwa, java, picard, pilon, samtools, ...}
# refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
    # # Absolute path to the H37RV genome
# scratchDir = "/n/scratch2/jy250/test100/"
# samFolder = "/n/scratch2/jy250/test100/dnaSamFiles/" 
# bamFolder = "/n/scratch2/jy250/test100/dnaBamFiles/" 
# cleanedFastqFolder = "/n/scratch2/jy250/test100/cleanedFastqs/"
# inputFolder = "/n/scratch2/jy250/test100/" # "/n/scratch2/jy250/22SampleInputs/"
# vcfFolder = "/n/scratch2/jy250/test100/vcfFiles/"
# slurmFolder = "/n/scratch2/jy250/test100/slurmLogs/"
# 

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
# lineageFolder = "/n/scratch2/jy250/lineageCalls"
# slurmFolder = "/n/scratch2/jy250/22SampleTest/slurmLogs/"


envName = "JerryEnv" # Environment with all of the packages pre-loaded
    # Packages included: {bwa, java, picard, pilon, samtools, ...}
refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
    # Absolute path to the H37RV genome
# scratchDir = "/n/scratch2/jy250/1000SampleTest"
samFolder = "/n/scratch2/jy250/1000SampleTest/samFiles/" 
bamFolder = "/n/scratch2/jy250/1000SampleTest/bamFiles/" 
cleanedFastqFolder = "/n/scratch2/jy250/1000SampleTest/cleanedFastqs/"
inputFolder = "/n/scratch2/jy250/eQTL_Analysis/19_7_31"  # "/n/scratch2/jy250/test100/" # 
vcfFolder = "/n/scratch2/jy250/1000SampleTest"
slurmFolder = "/n/scratch2/jy250/slurmLogs/"

countsFolder = "/n/scratch2/jy250/eQTL_Analysis/19_7_31/counts/"
lineageFolder = "/n/scratch2/jy250/lineageCalls"

experimentFolder = "/n/scratch2/jy250/eQTL_Analysis/19_7_31"
def genPairedFileName(fileName):
    """Creates the name of the paired end file, assuming that one is given a file that ends with "_1.fastq":
        XYZ_1.fastq => XYZ_2.fastq"""
    # TODO: Find a cleaner way to do this
    return fileName.replace("_1.fastq", "_2.fastq")
    
   
# def setUpSLURM(commands_list):
    # Nothing
    # Adds the heading of a slurm script, does nothing for now
    
def prepConda(commands_list, envName = envName):
    """Load the conda environment, which has already been created with a number of packages"""
    commands_list.append('module load conda2')
    commands_list.append('source deactivate') # Removes any pre-existing conda environments
    commands_list.append('source activate {eName}'.format(eName = envName))

def runQC(commands_list, fastq, fastqPaired = None, minQual = 10, sampleFolder = "./"):
    """ Quality control steps, including: {prinseq read trimming}"""
    # TODO: Add QC steps
    # TODO: Add Kraken
    return trimReads(commands_list, fastq, fastqPaired, minQual)

# def trimReads(commands_list, fastq, fastqPaired = None, minQual = 10, sampleFolder = "./"):
    # """Function to trim reads using fastp 
        # (https://github.com/farhat-lab/metatools_ncbi)
        
        # # Takes in a (pair of) FASTQ file(s) and appends to commands_list 
         # # the commands needed to trim each read of the files to a minimum phred 
         # # score of MINQUAL
    
        # # Outputs the name of the cleaned FASTQ file or files"""
    # leftMinQual = minQual
    # rightMinQual = minQual
    # sampleFolder = os.path.join(experimentFolder, genSampleID(fastq))
    # outputFile = os.path.join(sampleFolder, genBaseName(fastq) + ".cleaned") # Give it an absolute path # TODO: HARDCODED! BAD
        # #TODO: Hardcoded to write into the input folder => change later!!!
    # if fastqPaired is not None: # Paired end handling
        # trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                # -out_good {outputFileName} -fastq {fastq1} -fastq2 {fastqPaired}".format(leftMinQual = leftMinQual, 
                # rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq, fastqPaired = fastqPaired)
        # commands_list.append(trimCommand)
    # else:
        # trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                # -out_good {outputFileName} -fastq {fastq1}".format(leftMinQual = leftMinQual, 
                # rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq)
        # commands_list.append(trimCommand)    

        
def trimReads(commands_list, fastq, fastqPaired = None, minQual = 5):
    """Function to trim reads using prinseq
        
        Takes in a (pair of) FASTQ file(s) and appends to commands_list 
         the commands needed to trim each read of the files to a minimum phred 
         score of MINQUAL
         
        Outputs the name of the cleaned FASTQ file or files"""
         
    leftMinQual = minQual
    rightMinQual = minQual
    sampleFolder = os.path.join(experimentFolder, genSampleID(fastq))
    outputFile = os.path.join(sampleFolder, genBaseName(fastq) + ".cleaned") # Incomplete prefix of pprinseq to take in
    
    if fastqPaired is not None:
        trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                -out_good {outputFileName} -fastq {fastq1} -fastq2 {fastqPaired}".format(leftMinQual = leftMinQual, 
                rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq, fastqPaired = fastqPaired)
        commands_list.append(trimCommand)
        return genCleanedOutputName(outputFile, paired = False)
    else:
        trimCommand = "prinseq-lite.pl -trim_qual_left {leftMinQual} -trim_qual_right {rightMinQual} \
                -out_good {outputFileName} -fastq {fastq1}".format(leftMinQual = leftMinQual, 
                rightMinQual = rightMinQual, outputFileName = outputFile, fastq1= fastq)
        commands_list.append(trimCommand)
        return genCleanedOutputName(outputFile, paired = True)
    # head, tail = ntpath.split(outputFile)
    # fileName = tail or ntpath.basename(head)
    # return os.path.join(sampleFolder, fileName)
   
    
def genCleanedOutputName(outputFile, paired = False):
    """Generates the output filename that Prinseq would generate
        given an output file prefix 
        
        E.g. SRR1202546_1.cleaned =>
                    {SRR1202546_1.cleaned_1.fastq, SRR1202546_1.cleaned_2.fastq}"""
    if paired:
        return outputFile.replace("cleaned", "cleaned_1.fastq"), \
            outputFile.replace("cleaned", "cleaned_2.fastq")
    else:
        return outputFile.replace("cleaned", "cleaned_1.fastq"), ""
    # if fastq2 is not None:
        # return fastq1.replace(".fastq", ".cleaned.fastq"), fastq2.replace(".fastq", ".cleaned.fastq")
    # else:
        # return fastq1.replace(".fastq", ".cleaned.fastq"), ""

    
def convertToSAM(commands_list, refGenome, fastqFile, fastqPaired = "", numThreads = 1):
    """Converts the named FASTQ file to a SAM file
    
        Uses the bwa package, specifically the mem function"""
    #create SAM file
    samfile = genSamName(fastqFile) #generate the SAM file name
    # commands_list.append('bwa index {refGenome}'.format(refGenome = refGenome))
    commands_list.append("bwa mem -t {nThreads} -M {refGen} {fastq1} {fastq2}> {sam}".format(refGen = refGenome, fastq1 = fastqFile, fastq2 = fastqPaired, sam = samfile, nThreads = numThreads))
    return samfile

def convertToBAM(commands_list, samfile, numThreads = 1):
    """Converts the named SAM file to a BAM file
            Indexes and sorts as well
       Uses the samtools package
       
       RETURNS: the name of the new BAM file
    """
    bamfile = genBamName(samfile)
    commands_list.append('samtools faidx {refGenome}'.format(refGenome = refGenome))
    
    bamfile = genBamName(samfile)
    bamCommandsRaw = ["samtools view -S -b {samFile} > {bamFile}",
                    "samtools sort {samFile} -o {bamFile}",
                    "samtools index {bamFile}"] 
                    # With: 
                        #    {samFile} = location of samFile in current directory
                        #    {bamFile} = name of future bam file
    bamCommands = [line.format(samFile = samfile, bamFile = bamfile) for line in bamCommandsRaw]
        # Fill in the blanks
        
    commands_list.extend(bamCommands)
    return bamfile    

def convertToVCF(commands_list, bamFile, refGen = refGenome, outputDir = "./"):
        # TODO: can I add a default value that is a function of another argument?
    memoryReq = "200G"
    vcfName = genVCFName(bamFile)
    commands_list.append("pilon -Xmx{mem} --genome {ref} --bam {bam} --outdir {outDir} --vcf --output {outPrefix}".format(
                                                        mem = memoryReq, ref = refGen, bam = bamFile, outDir = outputDir, outPrefix =vcfName))
    return None

def submitSlurmScript(commands_list, outputName = None):
    """Submits the commands listed in commands_list to 
        the o2 cluster using Jerry's credentials and 
         default values of memory and time
      
    Default Memory: 60G
    Default Time: 12 hrs"""
    longString = ";".join(commands_list)
    print(longString.replace(";", "\n"))
    if outputName is not None:
        sCommand = 'sbatch -p short -c 1 -t 0-11:59 --mem=60G --mail-user=Jerry_Yang@hms.harvard.edu \
                --output {outputSlurm} --wrap="{commandString}"'.format(commandString = longString, outputSlurm = outputName)
    else: 
        sCommand = 'sbatch -p short -c 1 -t 0-11:59 --mem=60G --mail-user=Jerry_Yang@hms.harvard.edu \
                        --wrap="{0}"'.format(longString)
    os.system(sCommand)

def cleanBams(commands_list, bamfile):
    """Quality control: Uses PICARD to clean BAM of exact duplicates"""
    #   create BAM file with removed duplicates
    drbamfile = bamfile 
    # metricsFile = genMetricsName(drbamfile)
    # commands_list.append( "picard MarkDuplicates I={bam} O={dupRemBam} REMOVE_DUPLICATES=true M={metrics} ASSUME_SORT_ORDER=coordinate".format(bam = bamfile, dupRemBam = drbamfile, metrics = metricsFile))
    return drbamfile

def countReads(commands_list, bamFile,
        referenceGTF = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.gtf"):
    """Uses FeatureCounts to count the RNAseq reads from a BAM file """
    outputName = genReadcountsName(bamFile)
    countCommand = "featureCounts -t gene -g locus_tag -a {refGTF} -o {outputName} {bFile}".format(
                                               refGTF = referenceGTF, outputName = outputName, bFile = bamFile)
    commands_list.append(countCommand)
    
def callLineage(commands_list, vcfFile, snpScheme = "/n/data1/hms/dbmi/farhat/Jerry/fast-lineage-caller-vcf/snp_schemes/coll.tsv", outputFile = None):
    """Calls lineage using Luca's fast-lineage-caller-vcf script
        (https://github.com/farhat-lab/fast-lineage-caller-vcf)"""
    if outputFile == None:
        outputFile = genLinFileName(vcfFile)
    lineageCommand = "/n/data1/hms/dbmi/farhat/Jerry/fast-lineage-caller-vcf/bin/fast-lineage-caller-vcf.py {vcf} {scheme} > {out}".format(vcf = vcfFile, 
                                                                                                                            scheme = snpScheme, out = outputFile)    
    commands_list.append(lineageCommand)