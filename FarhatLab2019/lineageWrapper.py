from JerryPipeToolbox import callLineage, submitSlurmScript
import os 

# Script to run Luca Frechi's Lineage Calling script (fast-lineage-caller-vcf.py)
    # on all VCF files found in the specified input directory
    
inputFolder = "/n/scratch2/jy250/test100/"
for filename in os.listdir(inputFolder):
    if filename.endswith(".vcf"): 
        commands = []
        callLineage(commands, os.path.join(inputFolder, filename))
        print(commands[0])
        os.system(commands[0])