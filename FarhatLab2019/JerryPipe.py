import os
import pandas as pd
import numpy as np
# from slurmpy import Slurm
# import vcf
import shutil

###############################################################
########## Constants:##########################################
###############################################################
envName = "JerryEnv" # Environment with all of the packages pre-loaded
    # Packages included: {bwa, java, picard, pilon, samtools, ...}
refGenome = "/n/data1/hms/dbmi/farhat/Jerry/References/GCF_000195955.2_ASM19595v2_genomic.fasta"
    # Absolute path to the H37RV genome
scratchDir = "/n/scratch2/jy250/"
samFolder = "/n/scratch2/jy250/samFiles/" 
bamFolder = "/n/scratch2/jy250/bamFiles/" 
inputFolder = "/n/scratch2/jy250/fastqFiles/"
vcfFolder = "/n/scratch2/jy250/vcfFiles/"

# These 3 functions will take names of files, and generate the coresponding 
        # name of the new file
        # e.g. abc.fasta ==> abc.sam
def genSamName(fastq): # function to systematically generate absolute paths to 
                          # sam file names from FASTQ file names
    return os.path.join(samFolder, os.path.splitext(fastq)[0] + ".sam")
def genBamName(sam):
    return os.path.join(bamFolder, os.path.splitext(sam)[0] + ".bam")
def genVCFName(bam):
    return os.path.join(vcfFolder, os.path.splitext(bam)[0])
###############################################################
###############################################################

listOfFastqs = [] # ["SRR6176440_1.fastq", "SRR6176434_1.fastq"]
for filename in os.listdir(inputFolder):
    if filename.endswith(".fastq"): 
        listOfFastqs.append(filename)
        
tag = "tag"
output_dir = "/n/scratch2/jy250/output/"
scratch_dir = "/n/scratch2/jy250/"
   
# def setUpSLURM(commands_list):
    # Nothing
    # Adds the heading of a slurm script, does nothing for now
    
def prepConda(commands_list):
    """Load the conda environment, which has already been created with a number of packages"""
    commands_list.append('module load conda2')
    commands_list.append('source deactivate') # Removes any pre-existing conda environments
    commands_list.append('source activate {eName}'.format(eName = envName))

def runQC(commands_list):
    """ Quality control steps, including: {...}"""
     # TODO: Add QC steps
    # TODO: Add Kraken
    # TODO: Add read trimming

def convertToSAM(commands_list, refGenome, fastqFile, fastqPaired = ""):
    """Converts the named FASTQ file to a SAM file
    
        Uses the bwa package, specifically the mem function"""
    #create SAM file
    samfile = genSamName(fastqFile) #generate the SAM file name
    # commands_list.append('bwa index {refGenome}'.format(refGenome = refGenome))
    commands_list.append("bwa mem -p -M {refGen} {fastq1} {fastq2}> {sam}".format(refGen = refGenome, fastq1 = fastqFile, fastq2 = fastqPaired, sam = samfile))
    return samfile

def convertToBAM(commands_list, samfile):
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

def submitSlurmScript(commands_list):
    longString = ";".join(commands_list)
    sCommand = 'sbatch -p short -c 1 -t 0-11:59 --mem=50G --mail-user=Jerry_Yang@hms.harvard.edu --wrap="{0}"'.format(longString)
    os.system(sCommand)
    
## DEPRECATED: Decided not to remove duplicates from RNAseq data:
    ## https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/

# def cleanBams(commands_list, bamfile):
    # """Quality control: Uses PICARD to clean BAM of exact duplicates"""
    # #   create BAM file with removed duplicates
    # drbamfile = bamfile.replace(".bam", ".duprem.bam") # name for new, cleaned BAM file 
    # metricsFile = drbamfile[:-4]+'.metrics'
    # commands_list.append( "picard MarkDuplicates I={bam} O={dupRemBam} REMOVE_DUPLICATES=true M={metrics} ASSUME_SORT_ORDER=coordinate".format(bamfile, drbamfile, metrics = metricsFile)
    
    

    
def Launch_JerryPipe(fastqFile, fqf2 = None):
    '''
    This script is based on Max Marin's JankyPipe. Repurposed to work with Conda and RNAseq data.
    
    Converts a single (set of) FASTQ file into a VCF file
    
    "This script launches a job to call variants for the input fastq files against H37Rv
    using a number of packages. The important output (VCF, lineage info files, quality report)
    is stored in the output directory while the intermediary files (SAMs, trimmed fastqs, BAM, etc)
    are stored in a scratch directory."
    '''
    # Fundamental components:
        # - Generate a file (a SLURM script) with the following commands:
            # -- take a FASTQ file, and convert it to SAM via (bwa)
            # -- take the SAM file and convert it to BAM via (samtools)
            # -- *** take the SAM file and convert it to VCF via (Pilon) ***
            
        # - For all of the files, submit the job via SBATCH command
    
    
    # Store all commands in a list
    commands_list = []
    
    # setUpSLURM(commands_list) # Probably delete later, set up the requirements for the slurm script to run
    prepConda(commands_list)
    runQC(commands_list)
    
    # We just want to keep track of the names
    #TODO: Add functionality for paired end reads
    samFile = convertToSAM(commands_list, refGenome, fastqFile) # Adds commands to commandList to convert 
                                                      # the named fastq file to the named sam file
    bamFile = convertToBAM(commands_list, samFile)                 # Adds commands to commandList to convert 
                                                      # the named sam file to the named bam file
    vcfFile = convertToVCF(commands_list, bamFile, refGen = refGenome, outputDir = "/") 
    
    
    # "\n".join(commands_list)
    print(commands_list)
    submitSlurmScript(commands_list)
    # ######################################
    # ### QUALIMAP (quality of BAM file) ###
    # ######################################
    
    # #store quality report, pilon VCF & lineage call information all in Output directory
    # commands_list.append( 'cd ' + output_dir )
    # commands_list.append( 'mkdir QualiMap' ) #make a folder for pilon output in output directory
    # commands_list.append( 'unset DISPLAY' ) #unset JAVA virtual machine variable [http://qualimap.bioinfo.cipf.es/doc_html/faq.html]
    # commands_list.append( "/n/data1/hms/dbmi/farhat/bin/qualimap_v2.2.1/qualimap bamqc -bam {0} --outdir {1} --outfile {2}.pdf --outformat PDF".format(drbamfile, output_dir+'/QualiMap', tag+'_stats') )

    # ###################################
    # ### PILON (call variants) #########
    # ###################################
    
    # #store quality report, pilon VCF & lineage call information all in Output directory
    # commands_list.append( 'mkdir pilon' ) #make a folder for pilon output in output directory
    # out_pilon_dir = output_dir + '/pilon/' #variable for pilon output path

    # commands_list.append( 'java -Xmx32G -jar /n/data1/hms/dbmi/farhat/bin/pilon/pilon-1.22.jar --genome {0} --bam {1} --output {2} --outdir {3} --variant'.format(RefGen, drbamfile, tag, out_pilon_dir) )

    # #####################################
    # ### Luca's LINEAGE CALLING script ###
    # #####################################

    # #create directory 
    # commands_list.append( 'mkdir ' + scratch_dir + '/fast-lineage-caller/' )#make a folder for lineage call in output directory
    # commands_list.append( 'mkdir ' + output_dir + '/fast-lineage-caller/' )#make a folder for lineage call in scratch directory

    # #create VRT file
    # vrtfile = scratch_dir + '/fast-lineage-caller/{}.vrt'.format(tag)

    # commands_list.append( 'cd ' + scratch_dir + '/fast-lineage-caller' )#change directory to store output in scratch

    # #convert VCF to VRT
    # commands_list.append( '/n/data1/hms/dbmi/farhat/lfreschi/repos/vrt-tools/bin/vrtTools-vcf2vrt.py {0} {1} 1'.format(out_pilon_dir+tag+'.vcf', vrtfile) )

    # #call lineage with SNP database an VRT file
    # commands_list.append( 'cd ' + output_dir + '/fast-lineage-caller' )#change directory to store output in VCF output

    # commands_list.append( '/n/data1/hms/dbmi/farhat/lfreschi/repos/fast-lineage-caller/bin/FastLineageCaller-assign2lineage.py --lin_snps /home/rv76/Bio_Pipelines/fast-lineage-caller-master/example/db_snps.tsv ' + vrtfile + ' &> ' + 'lineage_call.txt' )

    # ###############################################################################################################
    # ######################################## SUBMIT as a job to O2 ################################################
    # ###############################################################################################################
    
    # #append all commands in a single string to be submitted as a job
    # JankyPipe_job = ''
    # for command_i in commands_list:
        # JankyPipe_job = JankyPipe_job + '\n' + command_i
        
        # #print(command_i)
        # #print(' ')
    
    
    # #directory where you want output + error files
    # # os.chdir(O2_SLURM_logs_dir)

    # job_name = tag

    # s = Slurm(job_name , {'partition':'short' , 'n':'1' , 't':'0-10:00:00' , 'mem-per-cpu':'48G' , 'mail-type':'FAIL' , 'mail-user':'mgmarin@g.harvard.edu'})

    # #submits the job
    # job_id = s.run(JankyPipe_job)

    # print(job_name  + ' : ' +  str(job_id))
for fastqFile in listOfFastqs:
    fileName = os.path.join(inputFolder, fastqFile)
    Launch_JerryPipe(fileName)
    
# Launch_JerryPipe("/n/scratch2/jy250/fastqFiles/SRR6176430_1.fastq")
