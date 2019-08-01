# Theoretically, this script should be able to be modified to run any command line exe file

import sys, os

# ONLY FOR COMPUTERS THAT ARE {1.) Linux, 2.) O2} 
cwd = os.getcwd() # Gets the current directory in O2 (under the assumption that my data folder will be relatively
                        # close to my current folder
                        

#outputFolder = "/n/data1/hms/dbmi/farhat/Jerry/inputs/"
# outputFolder = "/n/scratch2/jy250/sams/" # Send to much larger scratch folder
outputFolder = "" # Local to FarhatLab directory 
inputFolder = "/n/scratch2/jy250/"
condaEnv = "June27Env" # Name of Conda environment
slurmCommand = 'sbatch -p short -c 1 -t 0-11:59 --mem=50G --mail-user=Jerry_Yang@hms.harvard.edu --wrap="{0}"'
        # Base SLURM command to submit jobs
exe = "source activate {4}; bwa index {0}; bwa mem -M {1} {2} > {3}"
        # "samtools view -s -b {5}]
        # With: {0} = Reference genome
    #       {1} = Reference genome
    #       {2} = input .fastq file
    #       {3} = output .sam file
    #       {4} = conda environment name
    
bamToSam = " && samtools view -S -b {samFile} > {bamFile}; \
    samtools sort {samFile} -o {bamFile}; \
    samtools index {bamFile}"
        # With: 
        #    {samFile} = location of samFile in current directory
        #    {bamFile} = name of future bam file
    

refGenome = "../References/GCF_000195955.2_ASM19595v2_genomic.fasta"
        # bwa mem -M h37rv.fasta input/SRR6176429_1.fastq > testBWA.sam
        # idListLocation = "C:/Users/typew/Documents/GitHub/FarhatLab/3SamplePaperRunIDs.txt"
idListLocation = "/n/data1/hms/dbmi/farhat/Jerry/3SamplePaperRunIDs.txt" # Location of list of SRA IDs

# Now, want to loop thorugh every .fastq in the input directory and convert to BAM
for filename in os.listdir(inputFolder):
    if filename.endswith(".fastq"): 
        fileLocation = os.path.join(inputFolder, filename) # Gives the absolute path to the scratch2/ folder, where the FASTQ files are
       # Add the bamToSam commands
        outputLocation = os.path.join(outputFolder, filename + ".sam")
        bashCommand = exe.format(refGenome, refGenome, fileLocation, outputLocation, condaEnv)
        bamCommand = bamToSam.format(bamFile = filename + ".bam", samFile = filename + ".sam")

        sCommand = slurmCommand.format(bashCommand + bamCommand) # Bash command for SLURM to execute
        print(sCommand)
        os.system(sCommand)
        continue
    else:
        continue
  